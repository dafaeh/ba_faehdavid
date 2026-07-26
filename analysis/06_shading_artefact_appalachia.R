# 06_shading_artefact_appalachia.R
# This script works in the same way as the 05_shading_artefact_eel.R one but 
# for the appalachians focus region. 

# A shlightly different approach is used when filtering for high-quality pixels 
# only. Instead of using a Geology polygon to filter, the NLCD raster is used
# to filter for forest pixels.

# ---- Setup ------------------------------------------------------------------
library(terra)
library(tidyr)
library(dplyr)
library(ggplot2)
library(here)

source(here("R", "annotate_lm.R"))
source(here("R", "add_boxplot_n.R"))
source(here("R", "lanid_filter.R"))
source(here("R", "save_fig.R"))

set.seed(42)
n_sample <- 1000000

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
nlcd      <- stack_pre[["nlcd"]]    
lanid     <- stack_pre[["lanid"]]

# ---- Sample -----------------------------------------------------------------
stack <- c(cwd, elevation, twi, northness, slope)
names(stack) <- c("cwd_max", "elevation", "twi", "northness", "slope")

df_raw <- terra::spatSample(
  stack,
  size   = n_sample * 2,
  method = "regular",
  na.rm  = FALSE,
  xy     = TRUE,
  as.df  = TRUE
) |>
  tibble::as_tibble()

# LANID value per sampled pixel, used by the drop_irrigated() call below.
df_raw$lanid <- terra::extract(lanid, as.matrix(df_raw[, c("x", "y")]))$lanid

# Derive classification variables used in the plots below:
# aspect_class buckets northness into south-/east-west-/north-facing
# terciles (breaks at ±0.33). slope_class groups slope into gentle/moderate/steep bins
# elev_band splits elevation into deciles to relate slope to elevation.

# Minimum slope for a meaningful aspect. On near-flat terrain aspect is
# undefined and the algorithm assigns a random value, so northness
# carries no information there.
slope_min <- 2   # degrees

# Slope class labels are defined once here and reused by every scale that
# refers to them
slope_labels <- c(
  gentle   = paste0("gentle (", slope_min, " to 10 deg)"),
  moderate = "moderate (10 to 25 deg)",
  steep    = "steep (>25 deg)"
)

df <- df_raw |>
  drop_irrigated() |>
  tidyr::drop_na(cwd_max, elevation, twi, northness, slope) |>
  dplyr::filter(slope >= slope_min) |>
  mutate(
    aspect_class = cut(
      northness,
      breaks         = c(-1, -0.33, 0.33, 1),
      labels         = c("south-facing", "east/west", "north-facing"),
      include.lowest = TRUE
    ),
    slope_class = cut(
      slope,
      breaks         = c(slope_min, 10, 25, 90),
      labels         = unname(slope_labels),
      include.lowest = TRUE
    ),
    elev_band = cut(
      elevation,
      breaks         = quantile(elevation, probs = seq(0, 1, 0.1), na.rm = TRUE),
      labels         = paste0("Band ", 1:10),
      include.lowest = TRUE
    )
  )

# Sample size per aspect class, used to annotate Plots 2 and 3 below.
label_n_aspect_appalachia <- df |>
  dplyr::count(aspect_class) |>
  dplyr::mutate(line = paste0(aspect_class, ": n = ", format(n, big.mark = "'"))) |>
  dplyr::pull(line) |>
  paste(collapse = "\n")

# ---- Plot 1: CWDmax vs northness --------------------------------------------
# lm fit, matching the formula used by geom_smooth() below, so
# R²/n/slope can be read off with annotate_lm().
mod_northness_appalachia <- lm(cwd_max ~ northness, data = df)

plot_northness_appalachia <- ggplot(
  data = df,
  aes(x = northness, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "n (log10)") +
  geom_smooth(method = "lm", formula = y ~ x,
              se = FALSE, colour = "red", linewidth = 0.9) +
  scale_x_continuous(
    breaks       = c(-1, -0.5, 0, 0.5, 1),
    labels       = c("-1\nsouth", "-0.5", "0\neast/west", "0.5", "1\nnorth"),
    minor_breaks = NULL
  ) +
  labs(
    x = "Northness [cos(aspect)]",
    y = expression(CWD[max]~"[mm]")
  ) +
  theme_classic()

# Add R²/n/slope.
plot_northness_appalachia <- annotate_lm(
  plot_northness_appalachia, mod_northness_appalachia,
  slope_unit = "mm/unit northness"
)

save_fig(plot_northness_appalachia, "h2_appalachia_cwd_by_northness")

# ---- Plot 2: TWI vs CWDmax by aspect ----------------------------------------
plot_aspect_appalachia <- ggplot(
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

# Add per-aspect-group n label 
plot_aspect_appalachia <- plot_aspect_appalachia +
  annotate(
    "label",
    x = Inf, y = Inf, hjust = 1.05, vjust = 1.1,
    label = label_n_aspect_appalachia,
    size = 3.2, label.size = 0.3,
    label.padding = unit(0.4, "lines"),
    fill = scales::alpha("white", 0.9)
  )

save_fig(plot_aspect_appalachia, "h2_appalachia_cwd_aspect")

# ---- Plot 3: slope x aspect ------------------------------------------------
# Fixed y-axis limits, shared with Plot 3 in
# 06_shading_artefact_appalachia.R for cross-region comparability. 
shared_y_limits_slope_aspect <- c(137.7566, 826.7845)

plot_slope_aspect <- df |>
  group_by(aspect_class, slope_class) |>
  summarise(mean_cwd = mean(cwd_max, na.rm = TRUE), .groups = "drop") |>
  ggplot(
    aes(x = slope_class, y = mean_cwd, colour = aspect_class, group = aspect_class)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_colour_manual(values = aspect_colours, name = "Aspect") +
  scale_x_discrete(limits = rev) +
  coord_cartesian(ylim = shared_y_limits_slope_aspect) +
  labs(
    x       = "Slope class",
    y       = expression(CWD[max]~"[mm]")) +
  theme_classic() +
  theme(legend.position = "top")

# Per-aspect-group n label 
plot_slope_aspect <- plot_slope_aspect +
  annotate(
    "label",
    x = Inf, y = Inf, hjust = 1.05, vjust = 1.1,
    label = label_n_aspect_appalachia,
    size = 3.2, label.size = 0.3,
    label.padding = unit(0.4, "lines"),
    fill = scales::alpha("white", 0.9)
  )

save_fig(plot_slope_aspect, "h2_appalachia_slope_aspect")

# # ---- Plot 4: mean slope per elevation band ----------------------------------
# # To test whether CWDmax is systematically overestimated at higher elevations, 
# # which would confound H1.
# plot_slope_elev <- df |>
#   group_by(elev_band) |>
#   summarise(mean_slope = mean(slope, na.rm = TRUE), .groups = "drop") |>
#   ggplot(
#     aes(x = elev_band, y = mean_slope, group = 1)) +
#   geom_line(linewidth = 1, colour = "grey30") +
#   geom_point(size = 3, colour = "grey30") +
#   labs(
#     x       = "Elevation band (low to high)",
#     y       = "Mean slope [degrees]") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))
# 
# save_fig(plot_slope_elev, "h2_appalachia_slope_by_elev")
# 
# # ---- Plot 5: CWDmax by slope class ------------------------------------------
# plot_cwd_slope_class <- ggplot(
#   data = df, 
#   aes(x = slope_class, y = cwd_max, fill = slope_class)) +
#   geom_boxplot(outlier.alpha = 0.05, outlier.size = 0.5) +
#   # Names taken from slope_labels so the fill keys always match the factor
#   # levels created above, whatever slope_min is set to.
#   scale_fill_manual(
#     values = setNames(
#       c("#fee08b", "#fd8d3c", "#d7301f"),
#       unname(slope_labels)
#     )
#   ) +
#   labs(
#     x       = "Slope class",
#     y       = expression(CWD[max]~"[mm]")) +
#   theme_classic() +
#   theme(legend.position = "none")
# 
# # Add per-group n label.
# plot_cwd_slope_class <- plot_cwd_slope_class +
#   add_boxplot_n(df, "slope_class", "cwd_max")
# 
# save_fig(plot_cwd_slope_class, "h2_appalachia_cwd_by_slope")

# ---- Plot 6: TWI vs CWDmax, south-facing forest only ------------------------
# Isolate high-quality pixels: south-facing (least affected by the shading) 
# and forest only, to check whether a TWI-CWDmax signal emerges once aspect and land cover are held
# constant. Re-samples independently rather than reusing the shared df: over
# the full Appalachia stack as.data.frame() runs out of memory. slope is only
# carried along to reapply the slope_min filter.
stack_forest <- c(cwd, twi, northness, slope, nlcd, lanid)
names(stack_forest) <- c("cwd_max", "twi", "northness", "slope", "nlcd", "lanid")

df_forest_raw <- terra::spatSample(
  stack_forest,
  size   = n_sample * 2,
  method = "regular",
  na.rm  = FALSE,
  xy     = TRUE,
  as.df  = TRUE
) |>
  tibble::as_tibble()

# Remove irrigated pixels and filter for south facing forest pixels only
df_forest <- df_forest_raw |>
  drop_irrigated() |>
  dplyr::filter(nlcd %in% c(41, 42, 43)) |>   # deciduous/evergreen/mixed forest
  tidyr::drop_na(cwd_max, twi, northness, slope) |>
  dplyr::filter(slope >= slope_min) |>
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
mod_twi_forest_south_appalachia <- lm(cwd_max ~ twi, data = df_forest)

plot_twi_forest_south <- ggplot(
  data = df_forest,
  aes(x = twi, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "n (log10)") +
  geom_smooth(method = "lm", formula = y ~ x,
              se = FALSE, colour = "red", linewidth = 0.9) +
  labs(
    x = "TWI",
    y = expression(CWD[max]~"[mm]")) +
  theme_classic()

# Add R²/n/slope label
plot_twi_forest_south <- annotate_lm(
  plot_twi_forest_south, mod_twi_forest_south_appalachia,
  slope_unit = "mm/TWI unit"
)

save_fig(plot_twi_forest_south, "h2_appalachia_twi_forest_south")

# ---- Plot 7: TWI vs CWDmax, forest only, all aspects pooled -----------------
# Test whether the shading bias is also present, when all aspects are pooled.
# Same slope_min filter as everywhere else, so the pooled subset differs from
# Plot 6 only in the aspect filter.

stack_forest_pooled <- c(cwd, twi, northness, slope, nlcd, lanid)
names(stack_forest_pooled) <- c("cwd_max", "twi", "northness", "slope", "nlcd", "lanid")

df_forest_pooled_raw <- terra::spatSample(
  stack_forest_pooled,
  size   = n_sample * 2,
  method = "regular",
  na.rm  = FALSE,
  xy     = TRUE,
  as.df  = TRUE
) |>
  tibble::as_tibble()

# Remove irrigated pixels and filter for forest pixels only
df_forest_pooled <- df_forest_pooled_raw |>
  drop_irrigated() |>
  dplyr::filter(nlcd %in% c(41, 42, 43)) |> 
  tidyr::drop_na(cwd_max, twi, northness, slope) |>
  dplyr::filter(slope >= slope_min)

# Plot
mod_twi_forest_pooled_appalachia <- lm(cwd_max ~ twi, data = df_forest_pooled)

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

plot_twi_forest_pooled <- annotate_lm(
  plot_twi_forest_pooled, mod_twi_forest_pooled_appalachia,
  slope_unit = "mm/TWI unit"
)

save_fig(plot_twi_forest_pooled, "h2_appalachia_twi_forest_pooled")