# 05_shading_artefact_eel.R
# Provides evidence that a topographic shading artefact in DisALEXI prevents
# robust testing of H1 and H2 in the Eel River focus region.
#
# Restricted to eel_shading_subset, used to test sensitivity to shading
# effects on DisALEXI ET (see Risks and Contingency in the proposal).

# DisALEXI derives ET from land surface temperature (LST): low LST is
# interpreted as high ET via evaporative cooling. Topographic shading reduces
# LST independently of ET, causing ET and therefore CWDmax to be
# overestimated on shaded slopes. The artefact is expected to be strongest
# on steep north-facing slopes, where shading is strongest.

# ---- Setup ------------------------------------------------------------------
library(terra)
library(sf)
library(tidyr)
library(dplyr)
library(ggplot2)
library(here)

set.seed(42)
N_SAMPLE <- 1000000

# Colour scale shared across all plots
aspect_colours <- c(
  "south-facing" = "#d73027",
  "east/west"    = "#fee08b",
  "north-facing" = "#4575b4"
)

# Load data
cwd       <- terra::rast(here("data", "eel_9ref.tif"))[["cwd_max"]]

stack_pre <- terra::rast(here("data", "stack_eel.tif"))

# ---- Restrict to eel_shading_subset ------------------------------------------
aoi_wgs84 <- sf::st_polygon(list(matrix(
  c(-123.73946585911666, 39.301773660163114,
    -122.94845023411666, 39.301773660163114,
    -122.94845023411666, 41.2957488509508,
    -123.73946585911666, 41.2957488509508,
    -123.73946585911666, 39.301773660163114),
  ncol = 2, byrow = TRUE
))) |>
  sf::st_sfc(crs = "EPSG:4326")

aoi_5070 <- sf::st_transform(aoi_wgs84, "EPSG:5070")

cwd       <- terra::crop(cwd, terra::vect(aoi_5070))
stack_pre <- terra::crop(stack_pre, terra::vect(aoi_5070))

elevation <- stack_pre[["elevation"]]
twi       <- stack_pre[["twi"]]
northness <- stack_pre[["northness"]]
slope     <- stack_pre[["slope"]]
nlcd      <- stack_pre[["nlcd"]]
lanid     <- stack_pre[["lanid"]]

# ---- Sample -----------------------------------------------------------------
stack <- c(cwd, elevation, twi, northness, slope)
names(stack) <- c("cwd_max", "elevation", "twi", "northness", "slope")

df_raw <- terra::spatSample(
  stack,
  size   = N_SAMPLE * 2,
  method = "regular",
  na.rm  = FALSE,
  xy     = TRUE,
  as.df  = TRUE
) |>
  tibble::as_tibble()

# Exclude irrigated pixels using the LANID raster
df_raw$lanid <- terra::extract(lanid, as.matrix(df_raw[, c("x", "y")]))$lanid
# NLCD land cover, extracted here for the forest filter used in Plot 6.
df_raw$nlcd  <- terra::extract(nlcd, as.matrix(df_raw[, c("x", "y")]))$nlcd
n_irrigated  <- sum(df_raw$lanid == 1, na.rm = TRUE)

# Derive classification variables used in the plots below:
# aspect_class buckets northness into south-/east-west-/north-facing
# terciles (breaks at ±0.33) to test whether the shading artefact is
# directional. slope_class groups slope into gentle/moderate/steep bins
# to test whether the artefact scales with slope steepness. elev_band
# splits elevation into deciles to relate slope to elevation.
df <- df_raw |>
  dplyr::filter(is.na(lanid) | lanid != 1) |>
  tidyr::drop_na(cwd_max, elevation, twi, northness, slope) |>
  mutate(
    aspect_class = cut(
      northness,
      breaks         = c(-1, -0.33, 0.33, 1),
      labels         = c("south-facing", "east/west", "north-facing"),
      include.lowest = TRUE
    ),
    slope_class = cut(
      slope,
      breaks         = c(0, 10, 25, 90),
      labels         = c("gentle (<10 deg)", "moderate (10 to 25 deg)", "steep (>25 deg)"),
      include.lowest = TRUE
    ),
    elev_band = cut(
      elevation,
      breaks         = quantile(elevation, probs = seq(0, 1, 0.1), na.rm = TRUE),
      labels         = paste0("Band ", 1:10),
      include.lowest = TRUE
    )
  )

# ---- Plot 1: CWDmax vs northness --------------------------------------------
plot_northness_eel <- ggplot(
  data = df,
  aes(x = northness, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "n (log10)") +
  geom_smooth(method = "lm", formula = y ~ x,
              se = FALSE, colour = "red", linewidth = 0.9) +
  labs(
    x = "Northness (cos aspect) from south to north",
    y = expression(CWD[max]~"[mm]")) +
  theme_classic()

# Save plot
ggsave(
  filename = here("fig", "h2_eel_cwd_by_northness_subset.pdf"),
  plot     = plot_northness_eel,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_eel_cwd_by_northness_subset.png"),
  plot     = plot_northness_eel,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)

# ---- Plot 2: TWI vs CWDmax by aspect ----------------------------------------
plot_aspect_eel <- ggplot(
  data = df,
  aes(x = twi, y = cwd_max, colour = aspect_class)) +
  geom_point(alpha = 0.04, size = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  scale_colour_manual(values = aspect_colours, name = "Aspect") +
  labs(
    x = "TWI",
    y = expression(CWD[max]~"[mm]")) +
  theme_classic() +
  theme(legend.position = "top")

ggsave(
  filename = here("fig", "h2_eel_cwd_aspect_subset.pdf"),
  plot     = plot_aspect_eel,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_eel_cwd_aspect_subset.png"),
  plot     = plot_aspect_eel,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)

# ---- Plot 3: slope x aspect -------------------------------------------------
# Fixed y-axis limits for cross-region comparability with Plot 3 in
# 06_shading_artefact_appalachia.R. Eel River has the wider mean_cwd
# range, so it sets the bounds for both plots. Values from
# range(mean_cwd) below.
SHARED_Y_LIMITS_SLOPE_ASPECT <- c(137.7566, 810.8651)

plot_slope_aspect <- df |>
  group_by(aspect_class, slope_class) |>
  summarise(mean_cwd = mean(cwd_max, na.rm = TRUE), .groups = "drop") |>
  ggplot(
    aes(x = slope_class, y = mean_cwd, colour = aspect_class, group = aspect_class)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_colour_manual(values = aspect_colours, name = "Aspect") +
  scale_x_discrete(limits = rev) +
  coord_cartesian(ylim = SHARED_Y_LIMITS_SLOPE_ASPECT) +
  labs(
    x       = "Slope class",
    y       = expression(CWD[max]~"[mm]")) +
  theme_classic() +
  theme(legend.position = "top")

ggsave(
  filename = here("fig", "h2_eel_slope_aspect_subset.pdf"),
  plot     = plot_slope_aspect,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_eel_slope_aspect_subset.png"),
  plot     = plot_slope_aspect,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)

# ---- Plot 4: mean slope per elevation band ----------------------------------
# To test whether CWDmax is systematically overestimated at higher elevations, 
# which would confound H1.

plot_slope_elev <- df |>
  group_by(elev_band) |>
  summarise(mean_slope = mean(slope, na.rm = TRUE), .groups = "drop") |>
  ggplot(
    aes(x = elev_band, y = mean_slope, group = 1)) +
  geom_line(linewidth = 1, colour = "grey30") +
  geom_point(size = 3, colour = "grey30") +
  labs(
    x = "Elevation band (low to high)",
    y = "Mean slope (degrees)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))

ggsave(
  filename = here("fig", "h2_eel_slope_by_elev_subset.pdf"),
  plot     = plot_slope_elev,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_eel_slope_by_elev_subset.png"),
  plot     = plot_slope_elev,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)

# ---- Plot 5: CWDmax by slope class ------------------------------------------
plot_cwd_slope_class <- ggplot(
  data = df,
  aes(x = slope_class, y = cwd_max, fill = slope_class)) +
  geom_boxplot(outlier.alpha = 0.05, outlier.size = 0.5) +
  scale_fill_manual(
    values = c(
      "gentle (<10 deg)"        = "#fee08b",
      "moderate (10 to 25 deg)" = "#fd8d3c",
      "steep (>25 deg)"         = "#d7301f"
    )
  ) +
  labs(
    x = "Slope class",
    y = expression(CWD[max]~"[mm]")) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  filename = here("fig", "h2_eel_cwd_by_slope_subset.pdf"),
  plot     = plot_cwd_slope_class,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_eel_cwd_by_slope_subset.png"),
  plot     = plot_cwd_slope_class,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)

# ---- Plot 6: TWI vs CWDmax, south-facing forest only ------------------------
# Only high quality pixels are used for this plot. Restricting to south-facing
# pixels (the aspect class least affected by the shading bias) and to forest
# (NLCD 2019 codes 41/42/43) holds aspect and land cover constant.

df_south_forest <- df |>
  dplyr::filter(nlcd %in% c(41, 42, 43)) |>   # deciduous/evergreen/mixed forest
  dplyr::filter(aspect_class == "south-facing")

# Plot
plot_twi_forest_south <- ggplot(
  data = df_south_forest,
  aes(x = twi, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "n (log10)") +
  geom_smooth(method = "lm", formula = y ~ x,
              se = FALSE, colour = "red", linewidth = 0.9) +
  labs(
    x = "TWI",
    y = expression(CWD[max]~"[mm]")) +
  theme_classic()

# Save plot
ggsave(
  filename = here("fig", "h2_eel_twi_forest_south_subset.pdf"),
  plot     = plot_twi_forest_south,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_eel_twi_forest_south_subset.png"),
  plot     = plot_twi_forest_south,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)

# ---- Plot 7: TWI vs CWDmax, forest only, all aspects pooled -----------------
# Test whether the shading bias is also present, when all aspects are pooled.
# Reuses the shared df (already forest-filterable via nlcd), unlike the
# Appalachia counterpart which re-samples independently for the larger region.

df_forest_pooled <- df |>
  dplyr::filter(nlcd %in% c(41, 42, 43))   # deciduous/evergreen/mixed forest

# Plot
plot_twi_forest_pooled <- ggplot(
  data = df_forest_pooled,
  aes(x = twi, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "n (log10)") +
  geom_smooth(method = "lm", formula = y ~ x,
              se = FALSE, colour = "red", linewidth = 0.9) +
  labs(
    x = "TWI",
    y = expression(CWD[max]~"[mm]")) +
  theme_classic()

# Save plot
ggsave(
  filename = here("fig", "h2_eel_twi_forest_pooled_subset.pdf"),
  plot     = plot_twi_forest_pooled,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_eel_twi_forest_pooled_subset.png"),
  plot     = plot_twi_forest_pooled,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)