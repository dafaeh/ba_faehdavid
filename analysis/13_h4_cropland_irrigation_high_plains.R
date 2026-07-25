# 13_h4_cropland_irrigation_high_plains.R
# Tests H4 for the High Plains focus region: does CWDmax differ between
# irrigated and rainfed cropland?
# The classes are: Cropland irrigated and Cropland rainfed (NLCD 82 x LANID)
#
# LANID is the classification variable here, so no mask_irrigated() call:
# irrigated pixels are a class, not something to exclude. The ifel() chain
# below still drops pixels without LANID coverage (the condition evaluates to
# NA there), matching the project-wide rule that uncovered pixels are unusable.

# ---- Setup ------------------------------------------------------------------
library(terra)
library(dplyr)
library(ggplot2)
library(readr)
library(here)

# add_boxplot_n() adds per-group n to the boxplot below.
# save_fig() writes each figure to fig/ as PDF and PNG (see R/save_fig.R).
source(here("R", "add_boxplot_n.R"))
source(here("R", "save_fig.R"))

set.seed(42)

# Load CWD, NLCD and LANID

r     <- terra::rast(here("data", "high_plains_9ref.tif"))
r_cwd <- r[["cwd_max"]]

stack_pre <- terra::rast(here("data", "stack_high_plains.tif"))
nlcd      <- stack_pre[["nlcd"]]
lanid     <- stack_pre[["lanid"]]

# Combine NLCD and LANID into two classes Cropland irrigated and Cropland rainfed
cropland <- nlcd == 82

class_irrigation <- terra::ifel(
  cropland & lanid == 1, 1,
  terra::ifel(
    cropland & lanid == 0, 0,
    NA
  )
)
names(class_irrigation) <- "class"

# Stratified sample
# Irrigated and rainfed cropland end up comparably represented for the H4
# group comparison, despite different areal shares.
n_sample <- 100000

sample_pts <- terra::spatSample(
  class_irrigation,
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
      levels = c(1, 0),
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
    y = expression(CWD[max]~"[mm]")) +
  theme_classic() +
  theme(legend.position = "none")

# Add per-group n label.
plot_cropland_irrigation <- plot_cropland_irrigation +
  add_boxplot_n(df_filtered, "class", "cwd_max")

save_fig(plot_cropland_irrigation, "h4_cropland_irrigation_boxplot", width = 12, height = 10)