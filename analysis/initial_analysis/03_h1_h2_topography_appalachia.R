# 03_h1_h2_topography_appalachia.R
# Tests H1 (CWDmax decreases with elevation) and H2 (valley bottoms show
# higher CWDmax than ridges) for the Appalachian focus region.
#
# TWI is used as a proxy for topographic position: high TWI pixels
# correspond to convergent moisture-accumulating valley positions,
# low TWI pixels correspond to divergent ridge positions.
#
# H1: tested via regression of CWDmax against elevation and a mean
#     CWDmax per 100-m elevation band plot.
# H2: tested via regression of CWDmax against TWI.

# ---- Setup ------------------------------------------------------------------
library(terra)
library(dplyr)
library(ggplot2)
library(readr)
library(broom)
library(here)

set.seed(42)

path_cwd    <- here("data", "appalachia_9ref.tif")
path_elev   <- here("data", "elevation_appalachia.tif")
path_twi    <- here("data", "twi_appalachia.tif")
path_pixels <- here("data", "h1_h2_appalachia_pixels.csv")

n_sample    <- 100000
n_elev_band <- 10

# ---- Load CWD (needed before TWI crop) --------------------------------------
cwd <- terra::rast(path_cwd)[["cwd_max"]]

# ---- TWI: crop from CONUS file (cached) -------------------------------------
# The CONUS-wide TWI raster (EPSG:5072) is cropped to the study region and
# reprojected to EPSG:5070 to match the CWD grid. EPSG:5072 and EPSG:5070
# are both Conus Albers but use different datums (NAD83 NSRS2007 vs NAD83),
# so an explicit project() step is required.
# The result is cached to data/ so re-runs skip the expensive warp step.

if (!file.exists(path_twi)) {
  message("Preparing TWI for Appalachia ...")
  twi_conus <- terra::rast(
    here("data-raw", "CONUS_TWI_epsg5072_30m_unmasked.tif")
  )
  # Transform CWD extent to EPSG:5072 for efficient cropping of the CONUS file.
  ext_5072 <- terra::ext(
    terra::project(
      terra::as.polygons(terra::ext(cwd), crs = "EPSG:5070"),
      "EPSG:5072"
    )
  )
  twi_crop <- terra::crop(twi_conus, ext_5072)
  twi_proj <- terra::project(twi_crop, "EPSG:5070")
  twi      <- terra::resample(twi_proj, cwd, method = "bilinear")
  terra::writeRaster(twi, path_twi, overwrite = TRUE, datatype = "FLT4S")
  message("TWI saved to: ", path_twi)
} else {
  message("TWI cache found — skipping crop: ", path_twi)
  twi <- terra::rast(path_twi)
}

# ---- Load and align remaining layers ----------------------------------------
elev <- terra::resample(terra::rast(path_elev), cwd, method = "bilinear")
twi  <- terra::resample(twi, cwd, method = "bilinear")

stack        <- c(cwd, elev, twi)
names(stack) <- c("cwd_max", "elevation", "twi")

# ---- Stratified subsample ---------------------------------------------------
# Elevation bands are unequally distributed across the landscape. A simple
# random sample would undersample high-elevation pixels. Equal-sized strata
# based on elevation quantiles ensures all elevation ranges are represented.

df_full <- terra::spatSample(
  stack,
  size   = n_sample * 3,
  method = "regular",
  na.rm  = TRUE,
  as.df  = TRUE
)

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

readr::write_csv(df, path_pixels)
message("Pixel sample saved to: ", path_pixels)

# ---- H1: CWDmax ~ Elevation -------------------------------------------------
message("\n--- H1: CWDmax ~ elevation ---")

mod_h1 <- lm(cwd_max ~ elevation, data = df)
print(broom::tidy(mod_h1, conf.int = TRUE))
print(broom::glance(mod_h1)[c("r.squared", "adj.r.squared")])

# Mean CWDmax per 100-m elevation band reveals the structure of the
# relationship more clearly than a pixel-level scatter.
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

plot_elev_bands <- ggplot(df_h1_band, aes(x = mean_elev, y = mean_cwd)) +
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

ggsave(
  filename = here("fig", "h1_elev_bands_appalachia.pdf"),
  plot     = plot_elev_bands,
  width    = 16,
  height   = 12,
  units    = "cm"
)

# ---- H2: CWDmax ~ TWI -------------------------------------------------------
# High TWI corresponds to convergent valley positions with greater moisture
# accumulation; low TWI corresponds to divergent ridge positions.
# H2 predicts a positive relationship between TWI and CWDmax.

message("\n--- H2: CWDmax ~ TWI ---")

mod_h2 <- lm(cwd_max ~ twi, data = df)
print(broom::tidy(mod_h2, conf.int = TRUE))
print(broom::glance(mod_h2)[c("r.squared", "adj.r.squared")])

plot_cwd_twi <- ggplot(df, aes(x = twi, y = cwd_max)) +
  geom_point(alpha = 0.04, size = 0.3, colour = "grey40") +
  geom_smooth(method = "lm", colour = "firebrick", linewidth = 0.9) +
  labs(
    x     = "TWI",
    y     = "CWD_max [mm]",
    title = "CWD_max by TWI"
  ) +
  theme_classic()

ggsave(
  filename = here("fig", "h2_twi_appalachia.pdf"),
  plot     = plot_cwd_twi,
  width    = 16,
  height   = 12,
  units    = "cm"
)

message("\nAll outputs saved to fig/ and data/.")