# 02_h1_h2_topography_eel_river.R
# Tests H1 (CWDmax decreases with elevation) and H2 (valley bottoms show
# higher CWDmax than ridges at comparable altitude) for the Eel River focus
# region.
#
# TPI (Topographic Position Index) is used to classify ridge and valley pixels.
# TPI = elevation - mean elevation in a circular neighbourhood. Positive values
# indicate ridges, negative values indicate valleys. A 300 m radius is used to
# target hillslope-scale ridges and side valleys, which are the relevant scale
# for lateral subsurface flow (H2) and are only resolvable at 30 m resolution.
#
# H2 is tested by comparing CWDmax between valley and ridge pixels within
# elevation bands, so that the comparison is never confounded by elevation.


# ---- Setup ------------------------------------------------------------------

library(terra)
library(dplyr)
library(ggplot2)
library(readr)
library(broom)
library(here)

set.seed(42)

path_cwd    <- here("data", "eel_3ref.tif")
path_elev   <- here("data", "elevation_eel.tif")
path_tpi    <- here("data", "tpi_300_eel.tif")
path_pixels <- here("data", "h1_h2_eel_pixels.csv")

region      <- "eel"
tpi_radius  <- 300   # 300 meters should catch the hillslope scale
n_sample    <- 100000
n_elev_band <- 10


# ---- TPI computation (cached) -----------------------------------------------

# TPI is computed once and saved. On re-runs the cached file is loaded.

if (file.exists(path_tpi)) {
  message("TPI cache found, skipping computation: ", path_tpi)
} else {
  message("Computing TPI at ", tpi_radius, " m radius ...")
  w         <- terra::focalMat(terra::rast(path_elev), tpi_radius, type = "circle")
  elev_mean <- terra::focal(terra::rast(path_elev), w = w, fun = "mean", na.rm = TRUE)
  tpi       <- terra::rast(path_elev) - elev_mean
  names(tpi) <- "tpi_300"
  terra::writeRaster(tpi, path_tpi, overwrite = TRUE)
  message("TPI saved to: ", path_tpi)
}


# ---- Stack & stratified subsample -------------------------------------------

# Elevation bands are unequally distributed across the landscape. a simple
# random sample would undersample high-elevation pixels. Equal-sized strata
# based on elevation quantiles ensures all elevation ranges are represented.

cwd  <- terra::rast(path_cwd)[["cwd_max"]]
elev <- terra::rast(path_elev)
tpi  <- terra::rast(path_tpi)

elev <- terra::resample(elev, cwd, method = "bilinear")
tpi  <- terra::resample(tpi,  cwd, method = "bilinear")

stack        <- c(cwd, elev, tpi)
names(stack) <- c("cwd_max", "elevation", "tpi_300")

# Oversampling ensures enough valid pixels per elevation band after NA removal.
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

mod_h1 <- lm(cwd_max ~ elevation, data = df)
print(broom::glance(mod_h1))
print(broom::tidy(mod_h1, conf.int = TRUE))

# Mean CWDmax per 100-m elevation band reveals the structure of the
# relationship more clearly than a pixel-level scatter.
df_h1_band <- df |>
  dplyr::filter(!is.na(cwd_max)) |>
  dplyr::mutate(
    elev_band_100 = cut(
      elevation,
      breaks = seq(
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

p_h1 <- ggplot(df_h1_band, aes(x = mean_elev, y = mean_cwd)) +
  geom_ribbon(
    aes(ymin = mean_cwd - sd_cwd, ymax = mean_cwd + sd_cwd),
    fill = "grey40", alpha = 0.25
  ) +
  geom_line(linewidth = 1, colour = "grey20") +
  geom_point(size = 1.8, colour = "grey20") +
  labs(
    x       = "Elevation",
    y       = expression(mean~CWD[max]~"(mm)"),
    caption = "Mean \u00b1 1 SD per 100-m band"
  ) +
  theme_classic()

ggsave(
  filename = here("fig", "h1_elevation_bands_eel.pdf"),
  plot     = p_h1,
  width    = 16,
  height   = 12,
  units    = "cm"
)


# ---- H2: Valley vs. ridge within elevation bands ----------------------------

# Valley and ridge pixels are defined by the 25th and 75th TPI percentiles
# computed within each elevation band. This ensures that the comparison is
# always between pixels at similar altitudes. A ridge at 800 m is never
# compared to a valley at 200 m.

# Pixels in the middle 50% of the within-band TPI distribution are excluded
# to sharpen the contrast between convergent valley positions and divergent
# ridges.

df_h2 <- df |>
  dplyr::group_by(elev_band) |>
  dplyr::mutate(
    tpi_q25 = quantile(tpi_300, 0.25, na.rm = TRUE),
    tpi_q75 = quantile(tpi_300, 0.75, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(tpi_300 < tpi_q25 | tpi_300 > tpi_q75) |>
  dplyr::mutate(
    topo_class = dplyr::if_else(tpi_300 > tpi_q75, "Ridge", "Valley"),
    topo_class = factor(topo_class, levels = c("Valley", "Ridge"))
  )

# Elevation band labels: midpoint of each band rounded to nearest 50 m,
# formatted as a range string for the x-axis.
band_breaks <- quantile(
  df$elevation,
  probs          = seq(0, 1, 1 / n_elev_band),
  na.rm          = TRUE
)
band_labels <- sprintf(
  "%d\u2013%d m",
  round(band_breaks[-length(band_breaks)] / 50) * 50,
  round(band_breaks[-1]                  / 50) * 50
)

df_h2 <- df_h2 |>
  dplyr::mutate(
    elev_band_label = factor(
      band_labels[elev_band],
      levels = band_labels
    )
  )

p_h2 <- ggplot(df_h2, aes(x = elev_band_label, y = cwd_max, fill = topo_class)) +
  geom_boxplot(
    outlier.size  = 0.3,
    outlier.alpha = 0.15,
    linewidth     = 0.4
  ) +
  scale_fill_manual(
    values = c("Valley" = "grey40", "Ridge" = "grey80"),
    name   = NULL
  ) +
  labs(
    x = "Elevation band",
    y = expression(CWD[max]~"(mm)"),
  ) +
  theme_classic() +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 8),
    legend.position = "top"
  )

# ggsave(
#   filename = here("fig", "h2_valley_ridge_by_elev_band_eel.pdf"),
#   plot     = p_h2,
#   width    = 18,
#   height   = 12,
#   units    = "cm"
# )

message("\nAll outputs saved to fig/ and data/.")