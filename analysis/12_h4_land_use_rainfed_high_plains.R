# 12_h4_land_use_rainfed_high_plains.R
# Tests H4 for the High Plains focus region: does CWDmax differ across
# land-use types when irrigation is excluded?
# The classes are: Cropland rainfed, Grassland rainfed and Shrubland rainfed
# (NLCD 82 / 71 / 52 restricted to LANID = 0)


# ---- Setup ------------------------------------------------------------------
library(terra)
library(dplyr)
library(ggplot2)
library(here)

source(here("R", "add_boxplot_n.R"))
source(here("R", "lanid_filter.R"))
source(here("R", "save_fig.R"))

set.seed(42)

# Load CWD, NLCD and LANID
r <- terra::rast(here("data", "high_plains_9ref.tif"))

stack_pre <- terra::rast(here("data", "stack_high_plains.tif"))
nlcd      <- stack_pre[["nlcd"]]
lanid     <- stack_pre[["lanid"]]

# Reduce NLCD to three classes and keep only rainfed pixels.
# Pixels without LANID coverage are dropped too as their rainfed status is unknown.
nlcd_classes <- terra::classify(
  nlcd,
  rcl    = matrix(c(82, 1,
                    71, 2,
                    52, 3),
                  ncol = 2, byrow = TRUE),
  others = NA
)

class_rainfed <- mask_irrigated(nlcd_classes, lanid)
names(class_rainfed) <- "class"

# Stratified sample
# Cropland, Grassland, and Shrub/Scrub end up comparably represented for the
# group comparison
n_sample <- 100000

sample_pts <- terra::spatSample(
  class_rainfed,
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
      labels = c("Cropland rainfed", "Grassland rainfed", "Shrubland rainfed")
    )
  ) |>
  dplyr::filter(!is.na(cwd_max))


# Boxplot
plot_land_use_rainfed <- ggplot(
  data = df_filtered,
  aes(x = class, y = cwd_max)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, fill = "grey85") +
  labs(
    x = NULL,
    y = expression(CWD[max]~"[mm]")) +
  theme_classic() +
  theme(legend.position = "none")

# Add per-group n label.
plot_land_use_rainfed <- plot_land_use_rainfed +
  add_boxplot_n(df_filtered, "class", "cwd_max")

save_fig(plot_land_use_rainfed, "h4_2_land_use_rainfed", width = 12, height = 10)