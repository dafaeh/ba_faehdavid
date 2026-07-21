# 06_shading_artefact_appalachia.R
# This script works in the same way as the 05_shading_artefact_eel.R one but 
# for the appalachians focus region. It tests for the shading bias. 

# A shlightly different approach is used when filtering for high-quality pixels 
# only. Instead of using a Geology polygon to filter, the NLCD raster is used
# to filter for forest pixels.

# ---- Setup ------------------------------------------------------------------
library(terra)
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
cwd       <- terra::rast(here("data", "appalachia_9ref.tif"))[["cwd_max"]]

stack_pre <- terra::rast(here("data", "stack_appalachia.tif"))
elevation <- stack_pre[["elevation"]]
twi       <- stack_pre[["twi"]]
northness <- stack_pre[["northness"]]
slope     <- stack_pre[["slope"]]
nlcd      <- stack_pre[["nlcd"]]    # NLCD is needed later on to filter for high-quality pixels
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
n_irrigated  <- sum(df_raw$lanid == 1, na.rm = TRUE)

# Derive classification variables used in the plots below:
# aspect_class buckets northness into south-/east-west-/north-facing
# terciles (breaks at ±0.33) to test whether the shading artefact is
# directional. 
# slope_class groups slope into gentle/moderate/steep bins to test whether the 
# artefact scales with slope steepness. 
# elev_band splits elevation into deciles to relate slope to elevation.
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
plot_northness_appalachia <- ggplot(
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

ggsave(
  filename = here("fig", "h2_appalachia_cwd_by_northness.pdf"),
  plot     = plot_northness_appalachia,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_appalachia_cwd_by_northness.png"),
  plot     = plot_northness_appalachia,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)

# ---- Plot 2: TWI vs CWDmax by aspect ----------------------------------------
plot_aspect_appalachia <- ggplot(
  data = df, 
  aes(x = twi, y = cwd_max, colour = aspect_class)) +
  geom_point(alpha = 0.04, size = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  scale_colour_manual(values = aspect_colours, name = "Aspect") +
  labs(
    x     = "TWI",
    y       = expression(CWD[max]~"[mm]")) +
  theme_classic() +
  theme(legend.position = "top")

ggsave(
  filename = here("fig", "h2_appalachia_cwd_aspect.pdf"),
  plot     = plot_aspect_appalachia,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_appalachia_cwd_aspect.png"),
  plot     = plot_aspect_appalachia,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)

# ---- Plot 3: slope x aspect ------------------------------------------------
# Fixed y-axis limits for cross-region comparability with Plot 3 in
# 05_shading_artefact_eel.R. Bounds taken from Eel River's
# range(mean_cwd), which is wider than Appalachia's own range.
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
  filename = here("fig", "h2_appalachia_slope_aspect.pdf"),
  plot     = plot_slope_aspect,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_appalachia_slope_aspect.png"),
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
    x       = "Elevation band (low to high)",
    y       = "Mean slope [degrees]") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))

ggsave(
  filename = here("fig", "h2_appalachia_slope_by_elev.pdf"),
  plot     = plot_slope_elev,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_appalachia_slope_by_elev.png"),
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
    x       = "Slope class",
    y       = expression(CWD[max]~"[mm]")) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  filename = here("fig", "h2_appalachia_cwd_by_slope.pdf"),
  plot     = plot_cwd_slope_class,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_appalachia_cwd_by_slope.png"),
  plot     = plot_cwd_slope_class,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)

# ---- Plot 6: TWI vs CWDmax, south-facing forest only ------------------------
# Plot 2 showed south-facing pixels are the aspect class least affected by
# the shading artefact. This isolates that best-case subset further, to one
# land-cover class (deciduous/evergreen/mixed forest, NLCD 2019 codes
# 41/42/43), to check whether a TWI-CWDmax signal emerges once aspect and
# land cover are both held constant. Complete TWI range (no additional
# filtering).

# The Eel River counterpart (05_shading_artefact_eel.R, Plot 6) applies the
# same forest + south-facing filter, but reuses the shared df because its
# subset is small. Here the region is much larger: a full as.data.frame()
# over the Appalachia stack exhausts memory (std::bad_alloc), so this block
# re-samples independently (own stack incl. nlcd, own draw) instead. A regular
# sample avoids the memory blow-up and is dense enough for a hexbin. Sampling
# independently also leaves the shared df feeding Plots 1-5 untouched.

N_SAMPLE_FOREST <- 1000000

stack_forest <- c(cwd, twi, northness, nlcd, lanid)
names(stack_forest) <- c("cwd_max", "twi", "northness", "nlcd", "lanid")

df_forest_raw <- terra::spatSample(
  stack_forest,
  size   = N_SAMPLE_FOREST,
  method = "regular",
  na.rm  = FALSE,
  xy     = TRUE,
  as.df  = TRUE
) |>
  tibble::as_tibble()

# Remove irrigated pixels and filter for south facing forest pixels only
n_irrigated_forest <- sum(df_forest_raw$lanid == 1, na.rm = TRUE)

df_forest <- df_forest_raw |>
  dplyr::filter(is.na(lanid) | lanid != 1) |>
  dplyr::filter(nlcd %in% c(41, 42, 43)) |>   # deciduous/evergreen/mixed forest
  tidyr::drop_na(cwd_max, twi, northness) |>
  dplyr::mutate(
    aspect_class = cut(
      northness,
      breaks         = c(-1, -0.33, 0.33, 1),
      labels         = c("south-facing", "east/west", "north-facing"),
      include.lowest = TRUE
    )
  ) |>
  dplyr::filter(aspect_class == "south-facing")

# Plot
plot_twi_forest_south <- ggplot(
  data = df_forest,
  aes(x = twi, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "n (log10)") +
  geom_smooth(method = "lm", formula = y ~ x,
              se = FALSE, colour = "red", linewidth = 0.9) +
  labs(
    x     = "TWI",
    y       = expression(CWD[max]~"[mm]")) +
  theme_classic()

ggsave(
  filename = here("fig", "h2_appalachia_twi_forest_south.pdf"),
  plot     = plot_twi_forest_south,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_appalachia_twi_forest_south.png"),
  plot     = plot_twi_forest_south,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)

# ---- Plot 7: TWI vs CWDmax, forest only, all aspects pooled -----------------
# Test whether the shading bias is also present, when all aspects are pooled.

N_SAMPLE_FOREST_POOLED <- 1000000

stack_forest_pooled <- c(cwd, twi, northness, nlcd, lanid)
names(stack_forest_pooled) <- c("cwd_max", "twi", "northness", "nlcd", "lanid")

df_forest_pooled_raw <- terra::spatSample(
  stack_forest_pooled,
  size   = N_SAMPLE_FOREST_POOLED,
  method = "regular",
  na.rm  = FALSE,
  xy     = TRUE,
  as.df  = TRUE
) |>
  tibble::as_tibble()

# Remove irrigated pixels and filter for forest pixels only
n_irrigated_forest_pooled <- sum(df_forest_pooled_raw$lanid == 1, na.rm = TRUE)

df_forest_pooled <- df_forest_pooled_raw |>
  dplyr::filter(is.na(lanid) | lanid != 1) |>
  dplyr::filter(nlcd %in% c(41, 42, 43)) |> 
  tidyr::drop_na(cwd_max, twi, northness)

# Plot
plot_twi_forest_pooled <- ggplot(
  data = df_forest_pooled,
  aes(x = twi, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "n (log10)") +
  geom_smooth(method = "lm", formula = y ~ x,
              se = FALSE, colour = "red", linewidth = 0.9) +
  labs(
    x     = "TWI",
    y       = expression(CWD[max]~"[mm]")) +
  theme_classic()

ggsave(
  filename = here("fig", "h2_appalachia_twi_forest_pooled.pdf"),
  plot     = plot_twi_forest_pooled,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_appalachia_twi_forest_pooled.png"),
  plot     = plot_twi_forest_pooled,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)