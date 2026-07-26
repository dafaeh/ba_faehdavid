# plot_et_comparison.R
# Two side-by-side raster maps comparing OpenET Ensemble and DisALEXI
# annual-mean ET for the AOI exported by export_et_comparison.js.
# Shared colour scale (80 to 120 mm) taken from the etVis palette used
# for the GEE Code Editor inspection.


# ---- Setup --------------------------------------------------------------
library(terra)
library(tidyterra)
library(ggplot2)
library(ggspatial)
library(cowplot)
library(sf)
library(here)

# Load raster
r_et <- terra::rast(here("data-raw", "et_comparison_2021.tif"))

# Shared scale limits for optimal contrast
et_limits <- c(80, 120)

# Use a function to build the maps
make_et_map <- function(band_name, limits) {
  
  ggplot() +
    tidyterra::geom_spatraster(data = r_et, aes(fill = .data[[band_name]])) +
    scale_fill_viridis_c(
      name      = expression(ET~"(mm)"),
      limits    = limits,
      direction = -1,            # reversed: dark = high ET
      oob       = scales::squish,
      na.value  = NA
    ) +
    # Native CRS of the raster, no reprojection. datum = 4326 keeps the
    # graticule and axis labels in degrees.
    coord_sf(
      crs    = terra::crs(r_et),
      datum  = sf::st_crs(4326),
      expand = FALSE
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
        fill     = c("white", "grey40"),
        line_col = "grey20",
        text_col = "white"
      ),
      height = unit(0.8, "cm"),
      width  = unit(0.8, "cm")
    ) +

    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
    labs(x = "Longitude", y = "Latitude") +
    theme_classic(base_size = 13) +
    theme(
      axis.ticks        = element_line(colour = "grey30", linewidth = 0.4),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text         = element_text(size = 8),
      axis.text.x       = element_text(angle = 45, hjust = 1),
      axis.title        = element_text(size = 11)
    )
}

# Generate maps
map_disalexi <- make_et_map("et_disalexi", et_limits)
map_ensemble <- make_et_map("et_ensemble", et_limits)

# Display
map_disalexi
map_ensemble

# Save single panels
ggsave(here("fig", "map_et_disalexi.png"), plot = map_disalexi,
       width = 10, height = 8, units = "in", dpi = 600)
ggsave(here("fig", "map_et_ensemble.png"), plot = map_ensemble,
       width = 10, height = 8, units = "in", dpi = 600)

# Save combined 2-panel figure
fig_et_comparison <- plot_grid(
  map_disalexi, map_ensemble,
  labels = "auto", ncol = 2
)

ggsave(here("fig", "map_et_comparison.png"), plot = fig_et_comparison,
       width = 24, height = 10, units = "cm", dpi = 600)