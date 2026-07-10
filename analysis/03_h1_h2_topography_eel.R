# 03_h1_h2_topography_eel.R
# This Script tests H1 (CWDmax decreases with elevation) and H2 (valley bottoms show
# higher CWDmax than ridges) for the Eel River focus region.
#
# TWI is used as a proxy for topographic position: high TWI pixels
# correspond to convergent valley positions, where moisture accumulates while
# low TWI pixels correspond to divergent ridge positions.
#
# H1: tested via regression of CWDmax against elevation and a mean
#     CWDmax per 100m elevation band plot.
# H2: tested via regression of CWDmax against TWI.

# ---- Setup ------------------------------------------------------------------
library(terra)
library(dplyr)
library(ggplot2)
library(readr)
library(broom)
library(here)

set.seed(42) # Ensure reproducibility

# Load data
path_cwd    <- here("data", "eel_3ref.tif")
path_stack  <- here("data", "stack_eel.tif")
path_pixels <- here("data", "h1_h2_eel_pixels.csv")

# Define sample size and number of elevation bands
n_sample    <- 1000000
n_elev_band <- 10

# Load CWD raster 
cwd <- terra::rast(path_cwd)[["cwd_max"]]

# Load elevation and TWI from the preprocessed stack
stack_pre <- terra::rast(path_stack)
elev  <- stack_pre[["elevation"]]
twi   <- stack_pre[["twi"]]
lanid <- stack_pre[["lanid"]]   

stack        <- c(cwd, elev, twi)
names(stack) <- c("cwd_max", "elevation", "twi")

# ---- Stratified subsample ---------------------------------------------------
# Elevation bands are unequally distributed across the landscape. A simple
# random sample would undersample high-elevation pixels. Equal-sized strata
# based on elevation quantiles ensures all elevation ranges are represented.

df_full <- terra::spatSample(
  stack,
  size   = n_sample * 2,
  method = "regular",
  na.rm  = TRUE,
  xy     = TRUE,
  as.df  = TRUE
)

# Exclude irrigated pixels
df_full$lanid <- terra::extract(lanid, df_full[, c("x", "y")])$lanid
n_irrigated   <- sum(df_full$lanid == 1, na.rm = TRUE)
df_full       <- dplyr::filter(df_full, is.na(lanid) | lanid != 1)

# Quantile breaks are computed after the irrigation filter, so bands
# reflect the final population. Assumes each band retains at least
# n_per_band pixels post-filtering.
n_per_band <- n_sample / n_elev_band

df <- df_full |>
  dplyr::mutate(
    elev_band = cut(
      elevation,
      breaks         = quantile(elevation, probs = seq(0, 1, 1 / n_elev_band)),
      labels         = FALSE,
      include.lowest = TRUE
    )
  ) |>
  dplyr::group_by(elev_band) |>
  dplyr::slice_sample(n = n_per_band) |>
  dplyr::ungroup()

# ---- Test H1: CWDmax vs. Elevation ------------------------------------------
mod_h1 <- lm(cwd_max ~ elevation, data = df)
print(broom::tidy(mod_h1, conf.int = TRUE))
print(broom::glance(mod_h1)[c("r.squared", "adj.r.squared")])

df_h1_band <- df |>
  dplyr::filter(!is.na(cwd_max)) |>
  dplyr::mutate(
    elev_band_100 = cut(
      elevation,
      breaks         = seq(
        floor(min(elevation,   na.rm = TRUE)),
        ceiling(max(elevation, na.rm = TRUE)),
        by = 100
      ),
      include.lowest = TRUE
    )
  ) |>
  dplyr::group_by(elev_band_100) |>
  dplyr::summarise(
    mean_elev = mean(elevation, na.rm = TRUE),
    mean_cwd  = mean(cwd_max,   na.rm = TRUE),
    sd_cwd    = sd(cwd_max,     na.rm = TRUE),
    n         = sum(!is.na(cwd_max)),
    .groups   = "drop"
  ) |>
  dplyr::filter(n >= 50)

# Create a mean CWDmax per 100-m elevation band plot
plot_elev_bands <- ggplot(df_h1_band, aes(
  x = mean_elev, 
  y = mean_cwd)) +
  geom_ribbon(
    aes(ymin = mean_cwd - sd_cwd, ymax = mean_cwd + sd_cwd),
    fill = "grey40", alpha = 0.25
  ) +
  geom_line(linewidth = 1, colour = "grey20") +
  geom_point(size = 1.8, colour = "grey20") +
  labs(
    x       = "Elevation [m]",
    y       = "CWD_max [mm]",
    title   = "Mean CWD_max by Elevation",
    caption = "Mean and 1 SD per 100m band"
  ) +
  theme_classic()

# Save plot
ggsave(
  filename = here("fig", "h1_elev_bands_eel.pdf"),
  plot     = plot_elev_bands,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h1_elev_bands_eel.png"),
  plot     = plot_elev_bands,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)

# ---- Test H2: CWDmax vs. TWI ------------------------------------------------
mod_h2 <- lm(cwd_max ~ twi, data = df)
print(broom::tidy(mod_h2, conf.int = TRUE))
print(broom::glance(mod_h2)[c("r.squared", "adj.r.squared")])

# Plot CWDmax vs. TWI in scatterplot
plot_cwd_twi <- ggplot(df, aes(x = twi, y = cwd_max)) +
  geom_point(alpha = 0.04, size = 0.3, colour = "grey40") +
  geom_smooth(method = "lm", colour = "firebrick", linewidth = 0.9) +
  labs(
    x     = "TWI",
    y     = "CWD_max [mm]",
    title = "CWD_max by TWI"
  ) +
  theme_classic()

# Save plot
ggsave(
  filename = here("fig", "h2_twi_eel.pdf"),
  plot     = plot_cwd_twi,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_twi_eel.png"),
  plot     = plot_cwd_twi,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)