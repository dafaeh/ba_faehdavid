# 13_h4_cropland_irrigation_high_plains.R
# Tests H4 for the High Plains focus region: does CWDmax differ between
# irrigated and rainfed cropland?

# ---- Setup ------------------------------------------------------------------
library(terra)
library(dplyr)
library(ggplot2)
library(readr)
library(here)

set.seed(42)


# Load CWD raster
r     <- terra::rast(here("data", "high_plains_9ref.tif"))
r_cwd <- r[["cwd_max"]]


# Load NLCD and LANID
stack_pre     <- terra::rast(here("data", "stack_high_plains.tif"))
nlcd_aligned  <- stack_pre[["nlcd"]]
lanid_aligned <- stack_pre[["lanid"]]


# Combine NLCD and LANID into two classes Cropland irrigated and Cropland rainfed
cropland_mask <- nlcd_aligned == 82

class_img <- terra::ifel(
  cropland_mask & lanid_aligned == 1, 1,
  terra::ifel(
    cropland_mask & lanid_aligned == 0, 2,
    NA
  )
)
names(class_img) <- "class"


# Stratified sample
n_sample <- 500000

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


# ---- Quality filter ---------------------------------------------------------

df_filtered <- df_sample

# Eventually drop pixels with too many gap filled months, for example:

# df_filtered <- df_sample |>
#   dplyr::filter(gap_filled_months <= 6)

# Boxplot
plot_cropland_irrigation <- ggplot(
  data = df_filtered,
  aes(x = class, y = cwd_max)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, fill = "grey85") +
  labs(
    x = NULL,
    y = "CWD_max [mm]",
    title = "CWD_max Distribution by Irrigation Status within Cropland"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  here("fig", "h4_cropland_irrigation_boxplot.pdf"),
  plot   = plot_cropland_irrigation,
  width  = 12,
  height = 10,
  units  = "cm"
)

ggsave(
  here("fig", "h4_cropland_irrigation_boxplot.png"),
  plot   = plot_cropland_irrigation,
  width  = 12,
  height = 10,
  units  = "cm",
  dpi    = 600
)