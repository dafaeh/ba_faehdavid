# 05_shading_artefact_eel.R
# Provides evidence that a topographic shading artefact in DisALEXI prevents
# robust testing of H1 (CWDmax decreases with elevation) and H2 (valley
# bottoms show higher CWDmax than ridges) in the Eel River focus region.
#
# DisALEXI derives ET from land surface temperature (LST): low LST is
# interpreted as high ET via evaporative cooling. Topographic shading reduces
# LST independently of ET, causing ET and therefore CWDmax to be
# overestimated on shaded slopes. The artefact is expected to be strongest
# on steep north-facing slopes, where shading is most persistent.
#
# Argument structure:
#   1. Northness vs CWDmax: artefact is directional and significant
#   2. TWI vs CWDmax by aspect: artefact reverses the TWI signal, making
#      H2 untestable in the full dataset
#   3. Slope x aspect: artefact scales with slope steepness on north-facing
#      pixels, consistent with a geometric shading mechanism
#   4. Slope increases with elevation: same artefact confounds H1
#
# Nothing is saved. All output is printed to the console and plot window.

# ---- Setup ------------------------------------------------------------------
library(terra)
library(tidyr)
library(dplyr)
library(ggplot2)
library(here)

set.seed(42)
N_SAMPLE <- 150000

# Colour scale shared across all plots
aspect_colours <- c(
  "south-facing" = "#d73027",
  "east/west"    = "#fee08b",
  "north-facing" = "#4575b4"
)

# ---- Load and align rasters -------------------------------------------------
cwd       <- terra::rast(here("data", "eel_9ref.tif"))[["cwd_max"]]
elevation <- terra::resample(
  terra::rast(here("data", "elevation_eel.tif")), cwd, method = "bilinear"
)
twi <- terra::resample(
  terra::rast(here("data", "twi_eel.tif")), cwd, method = "bilinear"
)

# Terrain indices derived from the aligned elevation layer
northness <- cos(terra::terrain(elevation, "aspect", unit = "radians"))
names(northness) <- "northness"

slope <- terra::terrain(elevation, "slope", unit = "degrees")
names(slope) <- "slope"

# ---- Sample -----------------------------------------------------------------
stack <- c(cwd, elevation, twi, northness, slope)
names(stack) <- c("cwd_max", "elevation", "twi", "northness", "slope")

df <- terra::spatSample(
  stack,
  size   = N_SAMPLE,
  method = "regular",
  na.rm  = FALSE,
  as.df  = TRUE
) |>
  tibble::as_tibble() |>
  tidyr::drop_na() |>
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

message("Sample size after NA removal: ", nrow(df))

# ---- Plot 1: CWDmax vs northness --------------------------------------------
# The artefact is directional: CWDmax increases systematically with northness.
# R-squared is low because geology, vegetation and elevation also drive
# CWDmax variability. The coefficient (+96 mm per unit northness, p < 0.001)
# confirms a consistent bias regardless of low explained variance.

p1 <- ggplot(df, aes(x = northness, y = cwd_max)) +
  geom_point(alpha = 0.04, size = 0.3, colour = "grey40") +
  geom_smooth(method = "lm", colour = "firebrick", linewidth = 0.9) +
  coord_cartesian(ylim = c(0, 1000)) +
  labs(
    x       = "Northness (cos aspect)   south to north",
    y       = expression(CWD[max]~"(mm)"),
    title   = "CWDmax increases with northness (Eel River)",
    caption = "Coefficient: +96 mm per unit northness (p < 0.001)."
  ) +
  theme_classic()

print(p1)

# ---- Plot 2: TWI vs CWDmax by aspect ----------------------------------------
# The artefact does not merely shift CWDmax upward on north-facing slopes.
# It reverses the sign of the TWI signal: north-facing pixels show a negative
# TWI-CWDmax relationship, south-facing pixels show a positive one.
# The hydrological signal required to test H2 is therefore not recoverable
# from the full dataset.

p2 <- ggplot(df, aes(x = twi, y = cwd_max, colour = aspect_class)) +
  geom_point(alpha = 0.04, size = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  scale_colour_manual(values = aspect_colours, name = "Aspect") +
  coord_cartesian(ylim = c(0, 1000)) +
  labs(
    x     = "TWI",
    y     = expression(CWD[max]~"(mm)"),
    title = "TWI vs CWDmax by aspect (Eel River)"
  ) +
  theme_classic() +
  theme(legend.position = "top")

print(p2)

# ---- Plot 3: slope x aspect -------------------------------------------------
# If the artefact is caused by topographic shading, its magnitude should
# increase with slope on north-facing pixels and remain flat on south-facing
# pixels. This is the expected geometric signature of a shading mechanism.

p3 <- df |>
  group_by(aspect_class, slope_class) |>
  summarise(mean_cwd = mean(cwd_max, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(x = slope_class, y = mean_cwd,
             colour = aspect_class, group = aspect_class)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_colour_manual(values = aspect_colours, name = "Aspect") +
  labs(
    x       = "Slope class",
    y       = expression(mean~CWD[max]~"(mm)"),
    title   = "CWDmax increases with slope on north-facing pixels (Eel River)",
    caption = "Pattern is consistent with a geometric topographic shading artefact."
  ) +
  theme_classic() +
  theme(legend.position = "top")

print(p3)

# ---- Plot 4: mean slope per elevation band ----------------------------------
# Steep terrain is more frequent at higher elevations in the Eel River.
# Because the shading artefact is strongest on steep slopes, CWDmax is
# systematically overestimated at higher elevations, confounding H1 in the
# same direction as the expected signal.

p4 <- df |>
  group_by(elev_band) |>
  summarise(mean_slope = mean(slope, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(x = elev_band, y = mean_slope, group = 1)) +
  geom_line(linewidth = 1, colour = "grey30") +
  geom_point(size = 3, colour = "grey30") +
  labs(
    x       = "Elevation band (low to high)",
    y       = "Mean slope (degrees)",
    title   = "Slope increases with elevation (Eel River)",
    caption = "Steeper terrain at higher elevations amplifies the shading artefact, confounding H1."
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))

print(p4)

# ---- Regression summaries ---------------------------------------------------
message("\n--- lm: cwd_max ~ northness ---")
print(broom::tidy(lm(cwd_max ~ northness, data = df)))
print(broom::glance(lm(cwd_max ~ northness, data = df))[c("r.squared", "adj.r.squared")])

message("\n--- lm: cwd_max ~ twi ---")
print(broom::tidy(lm(cwd_max ~ twi, data = df)))

message("\n--- lm: cwd_max ~ twi + northness ---")
print(broom::tidy(lm(cwd_max ~ twi + northness, data = df)))
