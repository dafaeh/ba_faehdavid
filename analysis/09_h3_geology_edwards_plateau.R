# 09_h3_geology_edwards_plateau.R
# Tests H3 for the Edwards Plateau focus region.
# Compares CWDmax between carbonate and clastic geological substrates
# (SGMC field GENERALIZED_LITH).

# ---- Setup ------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(here)

set.seed(42)

# Load CWD and LANID rasters
r      <- terra::rast(here("data", "edwards_plateau_9ref.tif"))
r_cwd  <- r[["cwd_max"]]
r_gaps <- r[["gap_filled_months"]]

lanid <- terra::rast(here("data", "stack_edwards_plateau.tif"))[["lanid"]]

# Load and classify SGMC 
sgmc <- sf::st_read(
  here("data", "geology_edwards_plateau.gpkg"),
  quiet = TRUE
) |>
  dplyr::filter(GENERALIZED_LITH %in% c(
    "Sedimentary, carbonate",
    "Sedimentary, clastic"
  )) |>
  dplyr::mutate(
    geo_class = dplyr::case_when(
      GENERALIZED_LITH == "Sedimentary, carbonate" ~ "Carbonate",
      GENERALIZED_LITH == "Sedimentary, clastic"   ~ "Clastic"
    ),
    geo_class = factor(geo_class, levels = c("Clastic", "Carbonate"))
  )

# Rasterize the geology classes onto the CWD grid. This avoids loading all
# overlapping pixels into memory at once (which caused std::bad_alloc with
# terra::extract() on the large GENERALIZED_LITH polygons).
sgmc_vect <- terra::vect(sgmc)
sgmc_vect$class_int <- as.integer(sgmc$geo_class)

geo_rast <- terra::rasterize(sgmc_vect, r_cwd, field = "class_int")
names(geo_rast) <- "geo_class"

# Exclude irrigated pixels
geo_rast <- terra::mask(geo_rast, lanid, maskvalue = 1)

# Stratified sample 
# Per-class n_sample, not area-proportional: ensures Clastic and Carbonate
# are comparably represented for the H3 group comparison.
n_sample <- 500000

sample_pts <- terra::spatSample(
  geo_rast,
  size   = n_sample,
  method = "stratified",
  na.rm  = TRUE,
  xy     = TRUE,
  as.df  = TRUE
)

cwd_at_pts   <- terra::extract(r,     sample_pts[, c("x", "y")])
lanid_at_pts <- terra::extract(lanid, sample_pts[, c("x", "y")])

df_sample <- dplyr::tibble(
  geo_class         = sample_pts$geo_class,
  cwd_max           = cwd_at_pts$cwd_max,
  gap_filled_months = cwd_at_pts$gap_filled_months,
  lanid             = lanid_at_pts$lanid
) |>
  dplyr::mutate(
    geo_class = factor(
      geo_class,
      levels = c(1, 2),
      labels = levels(sgmc$geo_class)
    )
  ) |>
  dplyr::filter(!is.na(cwd_max))

# Boxplot 
p_box <- ggplot(
  df_sample,
  aes(x = geo_class, y = cwd_max, fill = geo_class)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
  scale_fill_manual(
    values = c("Clastic" = "#0072B2", "Carbonate" = "#E69F00")
  ) +
  labs(
    x = NULL,
    y = "CWD_max [mm]", 
    title = "CWD_max Distribution by Geological Substrate"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  here("fig", "h3_edwards_boxplot.pdf"),
  plot   = p_box,
  width  = 12,
  height = 10,
  units  = "cm"
)
ggsave(
  here("fig", "h3_edwards_boxplot.png"),
  plot   = p_box,
  width  = 12,
  height = 10,
  units  = "cm",
  dpi    = 600
)