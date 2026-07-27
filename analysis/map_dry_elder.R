# plot_dry_elder_sites.R
# Single satellite map showing both H3 validation field sites
# (Elder Creek, Dry Creek) with the catchment outlines.

# Polygons are delivered in NAD83 Contiguous USA Albers (parameter-identical to EPSG:5070) and are
# back-transformed to WGS84 only for display on the satellite tiles.

# ---- Setup --------------------------------------------------------------
library(maptiles)
library(terra)
library(tidyterra)
library(ggplot2)
library(ggspatial)
library(cowplot)
library(sf)
library(dplyr)
library(here)

# Panel size, same as the focus-region maps
panel_w_cm <- 16
panel_h_cm <- 12

# Catchment polygons
poly_elder <- sf::st_read(here("data-raw", "elder_creek_polygon.shp"), quiet = TRUE) |>
  dplyr::mutate(site = "Elder Creek")
poly_dry <- sf::st_read(here("data-raw", "dry_creek_polygon.shp"), quiet = TRUE) |>
  dplyr::mutate(site = "Dry Creek")

polys_5070 <- dplyr::bind_rows(poly_elder, poly_dry) |>
  dplyr::select(site)

polys_wgs <- sf::st_transform(polys_5070, crs = 4326)

# Label anchor points: st_point_on_surface() guarantees a point inside the
# polygon even for concave, hand-drawn shapes
labels_wgs <- sf::st_point_on_surface(polys_wgs)

# ---- Shared extent covering both sites --------------------------------------
# Bounding box of both catchments plus a margin (metres, EPSG:5070) so the
# polygons are not glued to the panel edge, then to WGS84 for the tiles.
margin_m <- 3000
box <- sf::st_bbox(polys_5070)
box["xmin"] <- box["xmin"] - margin_m
box["xmax"] <- box["xmax"] + margin_m
box["ymin"] <- box["ymin"] - margin_m
box["ymax"] <- box["ymax"] + margin_m

box_wgs <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(box), 4326))
xmin <- box_wgs[["xmin"]]; xmax <- box_wgs[["xmax"]]
ymin <- box_wgs[["ymin"]]; ymax <- box_wgs[["ymax"]]

# ---- Build the map ----------------------------------------------------------
ext   <- terra::ext(xmin, xmax, ymin, ymax)

# Vertical label nudge (map degrees), converted from ~1 cm on the saved
# figure using the y-extent and the ggsave() height 
fig_height_cm  <- panel_h_cm
deg_per_cm     <- (ymax - ymin) / fig_height_cm
nudge_y_labels <- 0.012 + 0.5 * deg_per_cm

tiles <- get_tiles(
  x        = ext,
  provider = "Esri.WorldImagery",
  zoom     = 13,
  crop     = TRUE,
  project  = FALSE
)

main_map <- ggplot() +
  geom_spatraster_rgb(data = tiles) +
  # Catchment outlines
  geom_sf(data = polys_wgs, fill = NA, colour = "red",
          linewidth = 0.4) +
  # Site labels next to the catchments
  geom_sf_text(data = labels_wgs, aes(label = site),
               colour = "white", size = 3, fontface = "bold",
               nudge_y = nudge_y_labels) +
# visible gap between the axis lines and the map
  coord_sf(
    xlim = c(xmin, xmax),
    ylim = c(ymin, ymax),
    crs  = 4326
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
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.ticks        = element_line(colour = "grey30", linewidth = 0.4),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text         = element_text(size = 8),
    # Rotated longitude labels take less horizontal space per label.
    # hjust = 1 keeps the label's end anchored under its tick mark.
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.title        = element_text(size = 11),
    plot.title        = element_text(size = 14)
  )

# Attribution as a cowplot label in the lower right corner, same position
# and styling as on the focus-region figures
map_sites <- ggdraw(main_map) +
  draw_label(
    "\u00a9 Esri World Imagery",
    x = 0.99, y = 0.01, hjust = 1, vjust = 0,
    size = 6, colour = "grey40"
  )


# Save
ggsave(here("fig", "map_dry_elder.png"), plot = map_sites,
       width = panel_w_cm, height = panel_h_cm, units = "cm", dpi = 600)
ggsave(here("fig", "map_dry_elder.pdf"), plot = map_sites,
       width = panel_w_cm, height = panel_h_cm, units = "cm", dpi = 600)