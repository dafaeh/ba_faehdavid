# 11_h4_land_use_high_plains.R
# Tests H4 for the High Plains focus region: does CWDmax differ across
# land-use types?
# The land-use types are: Cropland, Grassland and Shrubs
#
# No LANID filter at all here, unlike every other script: this is the
# all-land-use case over the full NLCD classes, and
# 12_h4_land_use_rainfed_high_plains.R is the rainfed counterpart.

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


# Load CWD and NLCD
r     <- terra::rast(here("data", "high_plains_9ref.tif"))
r_cwd <- r[["cwd_max"]]

stack_pre <- terra::rast(here("data", "stack_high_plains.tif"))
nlcd      <- stack_pre[["nlcd"]]


# Reduce NLCD to three classes Cropland, Grassland and Shrub
nlcd_aligned <- terra::classify(
  nlcd,
  rcl    = matrix(c(82, 1,
                    71, 2,
                    52, 3),
                  ncol = 2, byrow = TRUE),
  others = NA
)
names(nlcd_aligned) <- "class"


# Stratified sample
# Cropland, Grassland, and Shrub/Scrub end up comparably represented for the H4 
# group comparison, despite different areal shares.
n_sample <- 100000

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


# ---- Quality filter ---------------------------------------------------------

df_filtered <- df_sample

# Eventually drop pixels with too many gap filled months, for example:

# df_filtered <- df_sample |>
#   dplyr::filter(gap_filled_months <= 6)

# Boxplot
plot_land_use <- ggplot(
  data = df_filtered,
  aes(x = class, y = cwd_max)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, fill = "grey85") +
  labs(
    x = NULL,
    y = expression(CWD[max]~"[mm]")) +
  theme_classic() +
  theme(legend.position = "none")

# Add per-group n label.
plot_land_use <- plot_land_use + add_boxplot_n(df_filtered, "class", "cwd_max")

save_fig(plot_land_use, "h4_land_use_boxplot", width = 12, height = 10)