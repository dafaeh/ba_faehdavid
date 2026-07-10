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
# A land-cover-controlled comparison is provided in 09_h4_cropland_irrigation.R.
#
# Expected direction (H4): Irrigated > Non-irrigated.
# Rationale: irrigation supplements precipitation, allowing crops to
# sustain higher ET rates longer into the dry season, accumulating a
# larger CWD relative to non-irrigated land (proposal H4).
#
# Outputs:
#   - data/h4_test2_samples.csv
#   - fig/h4_irrigation_boxplot.pdf
#   - fig/h4_irrigation_boxplot.png


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

df_filtered <- df_sample

# Eventually drop pixels with too many gap-filled months, for example:

# df_filtered <- df_sample |>
#   dplyr::filter(gap_filled_months <= 6)


# ---- Boxplot ----------------------------------------------------------------

plot_irrigation <- ggplot(
  data = df_filtered,
  aes(x = class, y = cwd_max)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, fill = "grey85") +
  labs(
    x = NULL,
    y = "CWD_max [mm]",
    title = "CWD_max Distribution by Irrigation Status"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  here("fig", "h4_irrigation_boxplot.pdf"),
  plot   = plot_irrigation,
  width  = 12,
  height = 10,
  units  = "cm"
)

ggsave(
  here("fig", "h4_irrigation_boxplot.png"),
  plot   = plot_irrigation,
  width  = 12,
  height = 10,
  units  = "cm",
  dpi    = 600
)

message("\nFigures saved to fig/")