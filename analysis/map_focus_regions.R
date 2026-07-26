# plot_focus_regions.R
# Maps of the four BA focus regions and the Eel River subset

# Extent and CRS are read from the exported cwd_max rasters.
# Maps are drawn in EPSG:5070, datum = 4326 keeps
# the axis labels in degrees.


# ---- Setup --------------------------------------------------------------
library(maptiles)
library(terra)
library(tidyterra)
library(ggplot2)
library(ggspatial)   # annotation_scale(), annotation_north_arrow()
library(cowplot)
library(sf)
library(here)

# Define output sizes
panel_w_cm <- 16
panel_h_cm <- 12

# Use function to build maps
# Builds one focus-region map for every exported cwd_max
# raster.
make_region_map <- function(raster_path, title = NULL, zoom = 10,
                            show_title = TRUE, north_arrow = TRUE,
                            lon_break_width = 1, subset_bbox = NULL,
                            show_inset = TRUE, show_attribution = FALSE) {
  
  # Extent and CRS from the exported raster
  r       <- terra::rast(raster_path)
  crs_map <- sf::st_crs(terra::crs(r))
  
  # unname() is required
  e    <- as.vector(terra::ext(r))
  xmin <- unname(e["xmin"]); xmax <- unname(e["xmax"])
  ymin <- unname(e["ymin"]); ymax <- unname(e["ymax"])
  
  # Tiles are requested for a slightly larger area than the panel. The
  # basemap is warped from Web Mercator into EPSG:5070, and without the
  # margin that warp can leave a thin empty strip along the panel edge.
  buffer_m <- 0.02 * max(xmax - xmin, ymax - ymin)
  
  box_tiles <- sf::st_as_sfc(
    sf::st_bbox(
      c(xmin = xmin - buffer_m, ymin = ymin - buffer_m,
        xmax = xmax + buffer_m, ymax = ymax + buffer_m),
      crs = crs_map
    )
  )
  
  # project = TRUE reprojects the tiles into the CRS of x (EPSG:5070), so
  # the data layer is the one thing never resampled.
  tiles <- get_tiles(
    x        = box_tiles,
    provider = "Esri.WorldImagery",
    zoom     = zoom,
    crop     = TRUE,
    project  = TRUE
  )
  
  tiles <- terra::crop(tiles, terra::ext(xmin, xmax, ymin, ymax))
  
  main_map <- ggplot() +
    geom_spatraster_rgb(data = tiles)
  
  # Subset outline: drawn after the basemap so it sits on the satellite
  # image. 
  if (!is.null(subset_bbox)) {
    main_map <- main_map +
      annotate(
        "rect",
        xmin = subset_bbox[["xmin"]], xmax = subset_bbox[["xmax"]],
        ymin = subset_bbox[["ymin"]], ymax = subset_bbox[["ymax"]],
        fill = NA, colour = "red", linewidth = 0.9
      )
  }
  
  main_map <- main_map +
    # datum = 4326 keeps the axis labels in degrees. expand is left at its default (TRUE) so
    # ggplot pads a visible gap between the axis lines and the map
    coord_sf(
      crs   = crs_map,
      datum = sf::st_crs(4326),
      xlim  = c(xmin, xmax),
      ylim  = c(ymin, ymax)
    ) +
    
    annotation_scale(
      location   = "bl",
      width_hint = 0.15,
      style      = "bar",
      bar_cols   = c("black", "white"),
      height     = grid::unit(0.12, "cm"),
      text_cex   = 0.8
    ) +

    scale_x_continuous(breaks = scales::breaks_width(lon_break_width)) +
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
  
  if (show_title) {
    main_map <- main_map + labs(title = title)
  }
  
  if (north_arrow) {
    main_map <- main_map +
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
      )
  }
  
  # World map for inset. Skipped entirely when show_inset = FALSE, so
  # maps::map() isn't called for maps that don't need it.
  if (show_inset) {
    world <- sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))
    
    # The study box is drawn in the WGS84 inset, where the EPSG:5070
    # rectangle is a tilted quadrilateral. st_segmentize() adds vertices
    # along the edges before transforming, otherwise only the corners move
    # and the curved edges are drawn as straight lines.
    study_box <- sf::st_as_sfc(
      sf::st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
                  crs = crs_map)
    ) |>
      sf::st_segmentize(dfMaxLength = 5000) |>
      sf::st_transform(4326)
    
    inset_map <- ggplot() +
      geom_sf(data = world, fill = "grey30", colour = "grey50",
              linewidth = 0.1) +
      geom_sf(data = study_box, fill = "red", colour = "red", alpha = 0.4) +
      coord_sf(xlim = c(-130, -60), ylim = c(20, 55), expand = FALSE) +
      theme_void(base_size = 7) +
      theme(
        panel.background = element_rect(fill = "lightblue", colour = NA),
        panel.border     = element_rect(fill = NA, colour = "black",
                                        linewidth = 0.8)
      )
    
    final_map <- ggdraw(main_map) +
      draw_plot(inset_map, x = 0.71, y = 0.75, width = 0.20, height = 0.20)
  } else {
    final_map <- ggdraw(main_map)
  }
  
  if (show_attribution) {
    final_map <- final_map +
      draw_label(
        "\u00a9 Esri World Imagery",
        x = 0.99, y = 0.01, hjust = 1, vjust = 0,
        size = 6, colour = "grey40"
      )
  }
  
  final_map
}

# ---- Generate maps ----------------------------------------------------------
map_eel <- make_region_map(
  raster_path = here("data", "eel_9ref.tif"),
  title       = "Eel River",
  lon_break_width = 0.5
)
map_texas <- make_region_map(
  raster_path = here("data", "edwards_plateau_9ref.tif"),
  title       = "Edwards Plateau"
)
map_highplains <- make_region_map(
  raster_path = here("data", "high_plains_9ref.tif"),
  title       = "High Plains"
)
map_appalachians <- make_region_map(
  raster_path = here("data", "appalachia_9ref.tif"),
  title       = "Appalachians"
)

# Eel River sub-region
map_eel_subregion <- make_region_map(
  raster_path = here("data", "eel_9ref.tif"),
  lon_break_width = 0.5,
  subset_bbox = c(
    xmin = -2336730, xmax = -2229380,
    ymin =  2138000, ymax =  2366160
  ),
  show_inset       = FALSE,
  show_attribution = TRUE,
  show_title       = FALSE
)

# ---- Display ------------------------------------------------------------------
map_eel
map_texas
map_highplains
map_appalachians
map_eel_subregion

# ---- Save single panels --------------------------------------------------------
ggsave(here("fig", "map_eel_river.png"),    plot = map_eel,
       width = panel_w_cm, height = panel_h_cm, units = "cm", dpi = 600)
ggsave(here("fig", "map_texas.png"),        plot = map_texas,
       width = panel_w_cm, height = panel_h_cm, units = "cm", dpi = 600)
ggsave(here("fig", "map_highplains.png"),   plot = map_highplains,
       width = panel_w_cm, height = panel_h_cm, units = "cm", dpi = 600)
ggsave(here("fig", "map_appalachians.png"), plot = map_appalachians,
       width = panel_w_cm, height = panel_h_cm, units = "cm", dpi = 600)

ggsave(here("fig", "map_eel_subregion.png"), plot = map_eel_subregion,
       width = panel_w_cm, height = panel_h_cm, units = "cm", dpi = 600)

# ---- Save combined 4-panel figure ----------------------------------------------
fig_regions <- plot_grid(
  map_eel, map_texas, map_highplains, map_appalachians,
  labels = "auto", ncol = 2
)

fig_regions <- ggdraw(fig_regions) +
  draw_label(
    "\u00a9 Esri World Imagery",
    x = 0.99, y = 0.01, hjust = 1, vjust = 0,
    size = 6, colour = "grey40"
  )

# 2 x 2 panels at exactly the single-map size, so the scale bar, north
# arrow and axis text look the same here as on the single maps.
ggsave(here("fig", "map_focus_regions_combined.pdf"), plot = fig_regions,
       width = 2 * panel_w_cm, height = 2 * panel_h_cm, units = "cm",
       dpi = 600)

# PNG version of the same combined figure 
ggsave(here("fig", "map_focus_regions_combined.png"), plot = fig_regions,
       width = 2 * panel_w_cm, height = 2 * panel_h_cm, units = "cm",
       dpi = 600)