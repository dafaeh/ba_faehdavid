# 09_h4_cropland_irrigation_high_plains.R
# Tests H4 for the High Plains focus region: does CWDmax differ between
# irrigated and rainfed cropland? Combines NLCD (cropland mask) with
# LANID (irrigation status) to isolate the irrigation effect within a
# single land-cover class.
#
# Classification:
#   Cropland irrigated (class 1): NLCD 82 AND LANID = 1
#   Cropland rainfed   (class 2): NLCD 82 AND LANID = 0
#   All other pixels -> NA (non-cropland, irrigated pasture, etc.)
#
# By restricting both groups to NLCD class 82, differences in CWDmax
# are attributable to irrigation status rather than to contrasting
# vegetation types or land-management regimes.
#
# NLCD 2019 is streamed on demand via /vsicurl/. LANID 2017 must be
# present in data-raw/ (downloaded by 01_download_data.R).
#
# Expected direction (H4): Cropland irrigated > Cropland rainfed.
# Rationale: irrigation supplements precipitation, enabling irrigated
# crops to sustain higher ET further into the dry season and therefore
# accumulate a larger CWD than rainfed crops (proposal H4).
#
# Outputs:
#   - data/h4_test3_samples.csv
#   - fig/h4_cropland_irrigation_boxplot.pdf
#   - fig/h4_cropland_irrigation_density.pdf


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


# ---- Stream NLCD for AOI ----------------------------------------------------

# /vsicurl/ lets GDAL read the cloud-hosted GeoTIFF in place and only pull
# the tiles that overlap the AOI, avoiding a full CONUS download.
# Source: Dewitz (2023), https://doi.org/10.5066/P9JZ7AO3
nlcd_url  <- "/vsicurl/https://storage.googleapis.com/feddata-r/nlcd/2019_Land_Cover_L48.tif"
nlcd_full <- terra::rast(nlcd_url)

aoi_nlcd <- terra::project(
  terra::vect(terra::ext(r_cwd), crs = terra::crs(r_cwd)),
  terra::crs(nlcd_full)
)
nlcd <- terra::crop(nlcd_full, aoi_nlcd)


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


# ---- Align both rasters to the CWD pixel grid -------------------------------

# NLCD and LANID are both nominally EPSG:5070 at 30 m, but their pixel grids
# can be slightly offset. Nearest-neighbour resampling onto the CWD grid
# preserves categorical values and ensures pixel-by-pixel consistency.
nlcd_aligned  <- terra::resample(nlcd,  r_cwd, method = "near")
lanid_aligned <- terra::resample(lanid, r_cwd, method = "near")

# Minor CRS string differences between the two cause terra to warn.
terra::crs(lanid_aligned) <- terra::crs(nlcd_aligned)

cropland_mask <- nlcd_aligned == 82


# ---- Build cropland x irrigation class layer --------------------------------

# Class 1: cropland AND irrigated     (NLCD 82 AND LANID 1)
# Class 2: cropland AND non-irrigated (NLCD 82 AND LANID 0)
# All other pixels -> NA
cropland_mask <- nlcd_aligned == 82

class_img <- terra::ifel(
  cropland_mask & lanid_aligned == 1, 1,
  terra::ifel(
    cropland_mask & lanid_aligned == 0, 2,
    NA
  )
)
names(class_img) <- "class"


# ---- Stratified sample ------------------------------------------------------

# Draw 50 000 pixels per class. With millions of spatially autocorrelated
# pixels, a subsample is sufficient and keeps file sizes manageable.
n_sample <- 50000

sample_pts <- terra::spatSample(
  class_img,
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
      labels = c("Cropland irrigated", "Cropland rainfed")
    )
  ) |>
  dplyr::filter(!is.na(cwd_max))

message("Valid pixels by class:")
print(dplyr::count(df_sample, class))


# ---- Save subsampled data ---------------------------------------------------

readr::write_csv(df_sample, here("data", "h4_test3_samples.csv"))


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
  here("fig", "h4_cropland_irrigation_boxplot.pdf"),
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
  here("fig", "h4_cropland_irrigation_density.pdf"),
  plot   = p_dens,
  width  = 12,
  height = 8,
  units  = "cm"
)

message("\nFigures saved to fig/")