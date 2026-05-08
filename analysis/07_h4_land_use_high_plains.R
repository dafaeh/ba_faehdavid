# 07_h4_land_use_high_plains.R
# Tests H4 for the High Plains focus region: does CWDmax differ across
# land-use types?
#
# Land-cover classification (NLCD 2019, Dewitz 2023):
#   Cropland    (class 82): cultivated crops, the dominant agricultural use.
#   Grassland   (class 71): herbaceous grassland, the main natural cover.
#   Shrub/Scrub (class 52): shrubby vegetation, mainly in the drier west.
#
# Excluded classes and rationale:
#   Pasture/Hay (81) — too sparse (~1 % of AOI) and regionally specific
#                      (largely irrigated alfalfa); cannot be cleanly
#                      attributed to rainfed low-demand land use.
#   Forest, developed, water, wetlands — nearly absent in the High Plains
#                      and unrelated to H4.
#
# Expected direction (H4): Cropland > Grassland >= Shrub/Scrub.
# Rationale: cultivated crops have higher water demand than natural
# vegetation under otherwise comparable conditions (proposal H4).
#
# Outputs:
#   - data/h4_test1_samples.csv
#   - fig/h4_land_use_boxplot.pdf
#   - fig/h4_land_use_density.pdf


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

# CWD and NLCD use slightly different Albers CRS variants. Projecting the
# CWD extent into NLCD's CRS lets us crop NLCD on its native grid before
# any reprojection.
aoi_nlcd <- terra::project(
  terra::vect(terra::ext(r_cwd), crs = terra::crs(r_cwd)),
  terra::crs(nlcd_full)
)
nlcd <- terra::crop(nlcd_full, aoi_nlcd)


# ---- Reduce NLCD to three classes -------------------------------------------

# Of the 16 NLCD classes, only Cropland (82, ~48 % of AOI), Grassland (71,
# ~46 %) and Shrub/Scrub (52, ~5 %) cover meaningful fractions. All other
# classes are set to NA.
nlcd_class <- terra::classify(
  nlcd,
  rcl    = matrix(c(82, 1,
                    71, 2,
                    52, 3),
                  ncol = 2, byrow = TRUE),
  others = NA
)
names(nlcd_class) <- "class"

# Nearest-neighbour resampling is mandatory for categorical data.
nlcd_aligned <- terra::resample(nlcd_class, r_cwd, method = "near")


# ---- Stratified sample ------------------------------------------------------

# Draw 50 000 pixels per class. With millions of spatially autocorrelated
# pixels, a subsample is sufficient and keeps file sizes manageable.
n_sample <- 50000

sample_pts <- terra::spatSample(
  nlcd_aligned,
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
      levels = c(1, 2, 3),
      labels = c("Cropland", "Grassland", "Shrub/Scrub")
    )
  ) |>
  dplyr::filter(!is.na(cwd_max))

message("Valid pixels by class:")
print(dplyr::count(df_sample, class))


# ---- Save subsampled data ---------------------------------------------------

readr::write_csv(df_sample, here("data", "h4_test1_samples.csv"))


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
  here("fig", "h4_land_use_boxplot.pdf"),
  plot   = p_box,
  width  = 12,
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
  theme(legend.position = c(0.75, 0.80))

ggsave(
  here("fig", "h4_land_use_density.pdf"),
  plot   = p_dens,
  width  = 12,
  height = 8,
  units  = "cm"
)

message("\nFigures saved to fig/")