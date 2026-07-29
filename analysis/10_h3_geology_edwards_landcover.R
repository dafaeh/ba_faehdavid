# 10_h3_geology_edwards_landcover.R
# Extension of 09_h3_geology_edwards_plateau.R: checks whether
# the carbonate/clastic CWDmax contrast is hidden by the vegetation cover. 
# To test this, it is differentiated between the 3 dominant vegetation classes 
# on the Edwards Plateau: Grassland, Shrubland and Forest. 

# ---- Setup ------------------------------------------------------------------
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(here)

source(here("R", "lanid_filter.R"))
source(here("R", "save_fig.R"))

set.seed(42)

n_per_class <- 100000   # stratified sample per geology class (before veg split)


# Load CWD raster
r     <- terra::rast(here("data", "edwards_plateau_9ref.tif"))
r_cwd <- r[["cwd_max"]]


# Load, rasterize and classify geology 
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

sgmc_vect           <- terra::vect(sgmc)
sgmc_vect$class_int <- as.integer(sgmc$geo_class)
geo_rast            <- terra::rasterize(sgmc_vect, r_cwd, field = "class_int")
names(geo_rast)     <- "geo_class"


# Load NLCD
stack_pre    <- terra::rast(here("data", "stack_edwards_plateau.tif"))
nlcd_aligned <- stack_pre[["nlcd"]]

# Reclassify to three target classes; all others become NA.
# 1 = Grassland, 2 = Shrub, 3 = Forest
nlcd_reclass <- matrix(
  c(
    71, 1,          # Grassland/Herbaceous
    52, 2,          # Shrub/Scrub
    41, 3,          # Deciduous Forest
    42, 3,          # Evergreen Forest
    43, 3           # Mixed Forest
  ),
  ncol = 2, byrow = TRUE
)

veg_rast        <- terra::classify(nlcd_aligned, nlcd_reclass, others = NA)
names(veg_rast) <- "veg_class"


# Load LANID
lanid_aligned <- stack_pre[["lanid"]]


# Build combined class raster
# Mask geology to the three target vegetation classes and to pixels that
# LANID marks as non-irrigated.
combined <- terra::mask(geo_rast, veg_rast) # drop pixels outside the 3 veg classes
combined <- mask_irrigated(combined, lanid_aligned)
names(combined) <- "geo_class"


# ---- Stratified sample on geology class -------------------------------------
# n_per_class draws per geo_class stratum, not by area, so Clastic and
# Carbonate are comparably represented. veg_class is not itself stratified
sample_pts <- terra::spatSample(
  combined,
  size   = n_per_class,
  method = "stratified",
  na.rm  = TRUE,
  xy     = TRUE,
  as.df  = TRUE
)

ex_cwd <- terra::extract(r,        sample_pts[, c("x", "y")])
ex_veg <- terra::extract(veg_rast, sample_pts[, c("x", "y")])

veg_levels <- c("Grassland", "Shrub", "Forest")

df <- dplyr::tibble(
  geo_class         = sample_pts$geo_class,
  cwd_max           = ex_cwd$cwd_max,
  gap_filled_months = ex_cwd$gap_filled_months,
  veg_int           = ex_veg$veg_class
) |>
  dplyr::mutate(
    geo_class = factor(geo_class,
                       levels = c(1, 2),
                       labels = levels(sgmc$geo_class)),
    veg_class = factor(veg_int,
                       levels = seq_along(veg_levels),
                       labels = veg_levels)
  ) |>

  dplyr::filter(
    !is.na(cwd_max),
    !is.na(veg_class)
  )


# Faceted boxplot
# Same Wong-palette colours as 09_h3_geology_edwards_plateau.R so that
# geo_class is always encoded identically across figures.
geo_cols <- c("Clastic" = "#0072B2", "Carbonate" = "#E69F00")

p <- ggplot(df, aes(x = veg_class, y = cwd_max, fill = geo_class)) +
  geom_boxplot(
    outlier.size  = 0.3,
    outlier.alpha = 0.2,
    position      = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(values = geo_cols, name = NULL) +
  facet_wrap(~ veg_class, scales = "free_x", nrow = 1) +
  labs(
    x = NULL,
    y = expression(CWD[max]~"[mm]")
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    legend.position  = "bottom"
  )

# Per (veg_class x geo_class) n label.
df_n <- df |>
  dplyr::count(veg_class, geo_class)

p <- p +
  geom_text(
    data        = df_n,
    mapping     = aes(x = veg_class, y = max(df$cwd_max, na.rm = TRUE) * 1.05,
                      label = paste0("n = ", n), group = geo_class),
    position    = position_dodge(width = 0.8),
    inherit.aes = FALSE,
    size        = 2.6
  )

save_fig(p, "h3_edwards_landcover_boxplot", width = 16, height = 10)