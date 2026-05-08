# 08_h4_irrigation_high_plains.R
# Tests H4 for the High Plains focus region: does CWDmax differ between
# irrigated and non-irrigated land across all land-cover types?
#
# Irrigation classification (LANID 2017, Xie et al. 2021):
#   Irrigated     (LANID = 1): pixels mapped as irrigated in 2017.
#   Non-irrigated (LANID = 0): all remaining pixels (no irrigation signal).
#
# LANID 2017 is used as a proxy for irrigation status during the 2020-2022
# CWD analysis period. The dataset does not distinguish between land-cover
# types: irrigated pixels may include crops, hay, pasture, and golf courses.
# A land-cover-controlled comparison is provided in 07_h4_cropland_irrigation.R.
#
# Expected direction (H4): Irrigated > Non-irrigated.
# Rationale: irrigation supplements precipitation, allowing crops to
# sustain higher ET rates longer into the dry season, accumulating a
# larger CWD relative to non-irrigated land (proposal H4).
#
# Outputs:
#   - data/h4_test2_samples.csv
#   - fig/h4_irrigation_boxplot.pdf
#   - fig/h4_irrigation_density.pdf


# ---- Setup ------------------------------------------------------------------

library(terra)
library(dplyr)
library(ggplot2)
library(readr)
library(here)

set.seed(42)


# ---- Load raster ------------------------------------------------------------

stopifnot(
  "CWD raster not found. Place high_plains_9ref.tif in data/." =
    file.exists(here("data", "high_plains_9ref.tif"))
)

r     <- terra::rast(here("data", "high_plains_9ref.tif"))
r_cwd <- r[["cwd_max"]]


# ---- Load LANID -------------------------------------------------------------

# LANID is downloaded by 01_download_data.R and cached in data-raw/.
# Source: Xie et al. (2021), https://doi.org/10.5194/essd-13-5689-2021
stopifnot(
  "LANID 2017 not found. Run 01_download_data.R first." =
    file.exists(here("data-raw", "lanid2017.tif"))
)

lanid_full <- terra::rast(here("data-raw", "lanid2017.tif"))

aoi_lanid <- terra::project(
  terra::vect(terra::ext(r_cwd), crs = terra::crs(r_cwd)),
  terra::crs(lanid_full)
)
lanid <- terra::crop(lanid_full, aoi_lanid)


# ---- Build binary class layer -----------------------------------------------

# Retain only LANID values 0 and 1; set everything else to NA.
# Nearest-neighbour resampling preserves the binary classification.
lanid_aligned <- terra::resample(lanid, r_cwd, method = "near")

lanid_class <- terra::classify(
  lanid_aligned,
  rcl    = matrix(c(1, 1,
                    0, 2),
                  ncol = 2, byrow = TRUE),
  others = NA
)
names(lanid_class) <- "class"


# ---- Stratified sample ------------------------------------------------------

# Draw 50 000 pixels per class. With millions of spatially autocorrelated
# pixels, a subsample is sufficient and keeps file sizes manageable.
n_sample <- 50000

sample_pts <- terra::spatSample(
  lanid_class,
  size   = n_sample,
  method = "stratified",
  na.rm  = TRUE,
  xy     = TRUE,
  as.df  = TRUE
)

cwd_at_pts <- terra::extract(r, sample_pts[, c("x", "y")])

df_sample <- dplyr::tibble(
  x                 = sample_pts$x,
  y                 = sample_pts$y,
  class             = sample_pts$class,
  cwd_max           = cwd_at_pts$cwd_max,
  gap_filled_months = cwd_at_pts$gap_filled_months
) |>
  dplyr::mutate(
    class = factor(
      class,
      levels = c(1, 2),
      labels = c("Irrigated", "Non-irrigated")
    )
  ) |>
  dplyr::filter(!is.na(cwd_max))

message("Valid pixels by class:")
print(dplyr::count(df_sample, class))


# ---- Save subsampled data ---------------------------------------------------

readr::write_csv(df_sample, here("data", "h4_test2_samples.csv"))


# ---- Quality filter ---------------------------------------------------------

# Drop pixels where more than 6 of the 36 analysis months were gap-filled;
# these have a less reliable CWDmax.
df_filtered <- df_sample |>
  dplyr::filter(gap_filled_months <= 6)


# ---- Descriptive statistics -------------------------------------------------

df_desc <- df_filtered |>
  dplyr::group_by(class) |>
  dplyr::summarise(
    n          = dplyr::n(),
    median_mm  = median(cwd_max),
    mean_mm    = mean(cwd_max),
    sd_mm      = sd(cwd_max),
    q25_mm     = quantile(cwd_max, 0.25),
    q75_mm     = quantile(cwd_max, 0.75),
    gap_median = median(gap_filled_months, na.rm = TRUE),
    .groups    = "drop"
  )

message("\n=== Descriptive statistics ===")
print(df_desc, width = Inf)


# ---- Boxplot ----------------------------------------------------------------

p_box <- ggplot(
  df_filtered,
  aes(x = class, y = cwd_max, fill = class)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
  scale_fill_viridis_d(option = "D", end = 0.85) +
  labs(
    x = NULL,
    y = expression(paste(CWD[max], " (mm)"))
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  here("fig", "h4_irrigation_boxplot.pdf"),
  plot   = p_box,
  width  = 10,
  height = 10,
  units  = "cm"
)


# ---- Density plot -----------------------------------------------------------

p_dens <- ggplot(
  df_filtered,
  aes(x = cwd_max, fill = class, colour = class)) +
  geom_density(alpha = 0.4) +
  scale_fill_viridis_d(option = "D", end = 0.85, name = NULL) +
  scale_colour_viridis_d(option = "D", end = 0.85, name = NULL) +
  labs(
    x = expression(paste(CWD[max], " (mm)")),
    y = "Density"
  ) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.8))

ggsave(
  here("fig", "h4_irrigation_density.pdf"),
  plot   = p_dens,
  width  = 12,
  height = 8,
  units  = "cm"
)

message("\nFigures saved to fig/")