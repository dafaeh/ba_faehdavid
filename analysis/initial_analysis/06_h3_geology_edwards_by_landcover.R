# 06_h3_geology_edwards_landcover.R
# Exploratory extension of 06_h3_geology_edwards_plateau.R.
#
# Keeps the binary carbonate / clastic classification but adds two filters:
#   1. LANID mask: irrigated pixels are excluded, because supplemental
#      irrigation decouples CWDmax from natural water storage and would
#      confound the geology signal.
#   2. Vegetation class: pixels are grouped into three NLCD types that
#      dominate the Edwards Plateau under natural land cover:
#        Grassland  NLCD 71
#        Shrub      NLCD 52
#        Forest     NLCD 41, 42, 43
#      Cropland, developed, water, wetland are excluded. This lets us check
#      whether a geology effect is visible consistently across veg types,
#      or only within specific ones.
#
# Analysis: purely descriptive (n, median, IQR per geo x veg cell).
# A faceted boxplot shows the combined view.
#
# Inputs:
#   data/edwards_plateau_9ref.tif                CWD raster (cwd_max, gap_filled_months)
#   data-raw/USGS_StateGeologicMapCompilation_ver1.1.gdb   SGMC v1.1
#   data-raw/lanid2017.tif                       LANID irrigation layer
#   NLCD 2019 streamed via /vsicurl/
#
# Outputs:
#   data/h3_edwards_landcover_samples.csv
#   data/h3_edwards_landcover_desc.csv
#   fig/h3_edwards_landcover_boxplot.pdf


# ---- Setup ------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(here)

set.seed(42)

n_per_class <- 50000   # stratified sample per geology class (before veg split)
gap_max     <- 6       # max gap-filled months to keep a pixel


# ---- Load CWD raster --------------------------------------------------------

stopifnot(
  "CWD raster not found. Place edwards_plateau_9ref.tif in data/." =
    file.exists(here("data", "edwards_plateau_9ref.tif"))
)

r     <- terra::rast(here("data", "edwards_plateau_9ref.tif"))
r_cwd <- r[["cwd_max"]]


# ---- Load and rasterize SGMC (carbonate / clastic) --------------------------

# Identical to 06: read only carbonate and clastic polygons for the AOI.
aoi_wgs84 <- sf::st_as_sfc(sf::st_bbox(
  c(xmin = -100.27776, ymin = 28.78932,
    xmax =  -96.21282, ymax = 32.66018),
  crs = sf::st_crs(4326)
))

gdb_crs <- sf::st_crs(
  sf::st_read(
    here("data-raw", "USGS_StateGeologicMapCompilation_ver1.1.gdb"),
    query = "SELECT * FROM SGMC_Geology LIMIT 1",
    quiet = TRUE
  )
)

sgmc <- sf::st_read(
  dsn        = here("data-raw", "USGS_StateGeologicMapCompilation_ver1.1.gdb"),
  layer      = "SGMC_Geology",
  wkt_filter = sf::st_as_text(sf::st_transform(aoi_wgs84, gdb_crs)),
  quiet      = TRUE
) |>
  sf::st_transform(crs = 5070) |>
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


# ---- Stream and align NLCD --------------------------------------------------

nlcd_url  <- "/vsicurl/https://storage.googleapis.com/feddata-r/nlcd/2019_Land_Cover_L48.tif"
nlcd_full <- terra::rast(nlcd_url)

aoi_nlcd     <- terra::project(
  terra::vect(terra::ext(r_cwd), crs = terra::crs(r_cwd)),
  terra::crs(nlcd_full)
)
nlcd_aligned <- terra::resample(
  terra::crop(nlcd_full, aoi_nlcd),
  r_cwd,
  method = "near"
)

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


# ---- Load and align LANID ---------------------------------------------------

stopifnot(
  "LANID 2017 not found. Run 01_download_data.R first." =
    file.exists(here("data-raw", "lanid2017.tif"))
)

lanid_full <- terra::rast(here("data-raw", "lanid2017.tif"))

aoi_lanid      <- terra::project(
  terra::vect(terra::ext(r_cwd), crs = terra::crs(r_cwd)),
  terra::crs(lanid_full)
)
lanid_aligned  <- terra::resample(
  terra::crop(lanid_full, aoi_lanid),
  r_cwd,
  method = "near"
)
terra::crs(lanid_aligned) <- terra::crs(r_cwd)


# ---- Build combined class raster --------------------------------------------

# Mask: geology class is valid AND pixel is not irrigated AND veg class is valid.
# The LANID mask is applied by setting irrigated pixels to NA.
combined <- terra::ifel(
  lanid_aligned == 1,
  NA,              # exclude irrigated pixels
  geo_rast         # keep geo_class value for non-irrigated pixels
)
combined <- terra::mask(combined, veg_rast)   # drop pixels outside the 3 veg classes
names(combined) <- "geo_class"


# ---- Stratified sample on geology class -------------------------------------

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
    !is.na(veg_class),
    gap_filled_months <= gap_max
  )

message("Valid pixels per geo x veg cell:")
print(dplyr::count(df, geo_class, veg_class), n = Inf)

readr::write_csv(df, here("data", "h3_edwards_landcover_samples.csv"))


# ---- Descriptive statistics -------------------------------------------------

df_desc <- df |>
  dplyr::group_by(geo_class, veg_class) |>
  dplyr::summarise(
    n         = dplyr::n(),
    median_mm = median(cwd_max),
    mean_mm   = mean(cwd_max),
    sd_mm     = sd(cwd_max),
    q25_mm    = quantile(cwd_max, 0.25),
    q75_mm    = quantile(cwd_max, 0.75),
    .groups   = "drop"
  )

message("\n=== Descriptive statistics (mm) ===")
print(df_desc, width = Inf)

readr::write_csv(df_desc, here("data", "h3_edwards_landcover_desc.csv"))


# ---- Faceted boxplot --------------------------------------------------------

# Same Wong-palette colours as 06_h3_geology_edwards_plateau.R so that
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
    x       = NULL,
    y       = expression(CWD[max]~"(mm)"),
    title   = "CWD_max Distribution by Vegetation Class and Geological Substrate",
    caption = "Irrigated pixels are excluded."
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    legend.position  = "bottom"
  )

ggsave(
  here("fig", "h3_edwards_landcover_boxplot.pdf"),
  plot   = p,
  width  = 16,
  height = 10,
  units  = "cm"
)
ggsave(
  here("fig", "h3_edwards_landcover_boxplot.png"),
  plot   = p,
  width  = 16,
  height = 10,
  units  = "cm",
  dpi    = 600
)

message("\nOutputs written to data/ and fig/.")