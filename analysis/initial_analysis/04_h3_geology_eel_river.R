# 04_h3_geology_eel_river.R
# Tests H3 for the Eel River focus region: does CWDmax differ between
# geological units corresponding to the Coastal Belt (TK) and the
# Central Belt (KJf) of the Franciscan Complex?
#
# Geological classification (SGMC ORIG_LABEL):
#   TK  = "Tertiary-Cretaceous Coastal Belt Rocks"
#   KJf = "Franciscan Complex, unit 1 (Coast Ranges)"
#
# TK is explicitly named after the Coastal Belt. KJf corresponds
# geographically to the Central Belt in this region, as confirmed by
# visual inspection against Jayko et al. (1989) and Langenheim et al.
# (2013): KJf polygons lie consistently east of TK, separated by the
# Coastal Belt Thrust. Both units contain a single geological formation
# each (verified in sgmc_inventory_eel_river.csv).
#
# Expected direction (H3): Coastal Belt > Central Belt.
# Rationale: the Coastal Belt (TK) has deeper weathering profiles and
# greater subsurface water storage capacity than the Central Belt (KJf),
# enabling vegetation to sustain higher transpiration further into the
# dry season (Hahm et al. 2019, GRL).
#
# Outputs:
#   - data/h3_eel_river_pixels.csv
#   - fig/h3_eel_river_boxplot.pdf / .png
#   - fig/h3_eel_river_density.pdf / .png


# ---- Setup ------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(here)

set.seed(42)

# Wong-palette colours shared with 05_h3_cwd_dry_elder.R so that the
# Coastal Belt / Central Belt distinction is encoded identically across figures.
geo_cols <- c("Coastal Belt (TK)"  = "#009E73",
              "Central Belt (KJf)" = "#D55E00")


# ---- Load raster ------------------------------------------------------------

stopifnot(
  "CWD raster not found. Place eel_3ref.tif in data/." =
    file.exists(here("data", "eel_3ref.tif"))
)

r      <- terra::rast(here("data", "eel_3ref.tif"))
r_cwd  <- r[["cwd_max"]]
r_gaps <- r[["gap_filled_months"]]


# ---- Load and classify SGMC -------------------------------------------------

aoi_wgs84 <- sf::st_as_sfc(sf::st_bbox(
  c(xmin = -123.74888, ymin = 39.29197,
    xmax = -121.09019, ymax = 41.29345),
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
  dplyr::filter(ORIG_LABEL %in% c("TK", "KJf")) |>
  dplyr::mutate(
    geo_class = dplyr::case_when(
      ORIG_LABEL == "TK"  ~ "Coastal Belt (TK)",
      ORIG_LABEL == "KJf" ~ "Central Belt (KJf)"
    ),
    geo_class = factor(geo_class,
                       levels = c("Coastal Belt (TK)", "Central Belt (KJf)"))
  )


# ---- Extract pixels ---------------------------------------------------------

# Extract both layers to allow quality filtering via gap_filled_months
# downstream if needed.
extracted_cwd  <- terra::extract(r_cwd,  terra::vect(sgmc), fun = NULL)
extracted_gaps <- terra::extract(r_gaps, terra::vect(sgmc), fun = NULL)

df_pixels <- dplyr::tibble(
  ID                = extracted_cwd$ID,
  cwd_max           = extracted_cwd$cwd_max,
  gap_filled_months = extracted_gaps$gap_filled_months,
  geo_class         = sgmc$geo_class[extracted_cwd$ID],
  orig_label        = sgmc$ORIG_LABEL[extracted_cwd$ID]
) |>
  dplyr::filter(!is.na(cwd_max))

message("Valid pixels by class:")
print(dplyr::count(df_pixels, geo_class))


# ---- Subsample --------------------------------------------------------------

# With millions of spatially autocorrelated pixels, a subsample is
# sufficient and keeps file sizes manageable.
n_sample <- 50000

df_sample <- df_pixels |>
  dplyr::group_by(geo_class) |>
  dplyr::slice_sample(n = n_sample) |>
  dplyr::ungroup()


# ---- Save subsampled data ---------------------------------------------------

readr::write_csv(df_sample, here("data", "h3_eel_river_pixels.csv"))


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
  scale_fill_manual(values = geo_cols) +
  labs(
    x     = NULL,
    y     = expression(CWD[max]~"(mm)"),
    title = "CWD_max Distribution by Geological Formation"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  here("fig", "h3_eel_river_boxplot.pdf"),
  plot   = p_box,
  width  = 12,
  height = 10,
  units  = "cm"
)
ggsave(
  here("fig", "h3_eel_river_boxplot.png"),
  plot   = p_box,
  width  = 12,
  height = 10,
  units  = "cm",
  dpi    = 600
)

message("\nFigures saved to fig/")