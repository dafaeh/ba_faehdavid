# 06_h3_geology_edwards_plateau.R
# Tests H3 for the Edwards Plateau focus region.
# Compares CWDmax between carbonate and clastic geological substrates
# (SGMC field GENERALIZED_LITH).
# H3 predicts higher CWDmax in clastic substrates due to deeper soils
# and greater water retention relative to karstified carbonates.
#
# Outputs:
#   - data/h3_edwards_pixels.csv
#   - fig/h3_edwards_boxplot.pdf
#   - fig/h3_edwards_density.pdf


# ---- Setup ------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(here)

set.seed(42)


# ---- Load raster ------------------------------------------------------------

stopifnot(
  "CWD raster not found. Place edwards_plateau_9ref.tif in data/." =
    file.exists(here("data", "edwards_plateau_9ref.tif"))
)

r      <- terra::rast(here("data", "edwards_plateau_9ref.tif"))
r_cwd  <- r[["cwd_max"]]
r_gaps <- r[["gap_filled_months"]]


# ---- Load and classify SGMC -------------------------------------------------

stopifnot(
  "SGMC not found. Run 01_download_data.R first." =
    file.exists(here("data-raw", "USGS_StateGeologicMapCompilation_ver1.1.gdb"))
)

aoi_wgs84 <- sf::st_as_sfc(sf::st_bbox(
  c(xmin = -100.27776, ymin = 28.78932,
    xmax =  -96.21282, ymax = 32.66018),
  crs = sf::st_crs(4326)
))

# Read one row first to get the native GDB CRS, then use it for the
# spatial filter — avoids a full table scan on the large geodatabase.
gdb_crs <- sf::st_crs(
  sf::st_read(
    here("data-raw", "USGS_StateGeologicMapCompilation_ver1.1.gdb"),
    query = "SELECT * FROM SGMC_Geology LIMIT 1",
    quiet = TRUE
  )
)

aoi_gdb <- sf::st_transform(aoi_wgs84, crs = gdb_crs)

sgmc <- sf::st_read(
  dsn        = here("data-raw", "USGS_StateGeologicMapCompilation_ver1.1.gdb"),
  layer      = "SGMC_Geology",
  wkt_filter = sf::st_as_text(aoi_gdb),
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


# ---- Rasterize SGMC ---------------------------------------------------------

# Rasterize the geology classes onto the CWD grid. This avoids loading all
# overlapping pixels into memory at once (which caused std::bad_alloc with
# terra::extract() on the large GENERALIZED_LITH polygons).
# Nearest-neighbour is mandatory for categorical data.
sgmc_vect <- terra::vect(sgmc)
sgmc_vect$class_int <- as.integer(sgmc$geo_class)

geo_rast <- terra::rasterize(sgmc_vect, r_cwd, field = "class_int")
names(geo_rast) <- "geo_class"


# ---- Stratified sample ------------------------------------------------------

n_sample <- 50000

sample_pts <- terra::spatSample(
  geo_rast,
  size   = n_sample,
  method = "stratified",
  na.rm  = TRUE,
  xy     = TRUE,
  as.df  = TRUE
)

cwd_at_pts <- terra::extract(r, sample_pts[, c("x", "y")])

df_sample <- dplyr::tibble(
  geo_class         = sample_pts$geo_class,
  cwd_max           = cwd_at_pts$cwd_max,
  gap_filled_months = cwd_at_pts$gap_filled_months
) |>
  dplyr::mutate(
    geo_class = factor(
      geo_class,
      levels = c(1, 2),
      labels = levels(sgmc$geo_class)
    )
  ) |>
  dplyr::filter(!is.na(cwd_max))

message("Valid pixels by class:")
print(dplyr::count(df_sample, geo_class))


# ---- Save subsampled data ---------------------------------------------------

readr::write_csv(df_sample, here("data", "h3_edwards_pixels.csv"))


# ---- Descriptive statistics -------------------------------------------------

df_desc <- df_sample |>
  dplyr::group_by(geo_class) |>
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

message("\nFigures saved to fig/")