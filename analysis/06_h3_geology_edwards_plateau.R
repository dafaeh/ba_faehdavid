# 06_h3_geology_edwards_plateau.R
# Tests H3 for the Edwards Plateau focus region: does CWDmax differ
# between carbonate and non-carbonate geological substrates?
#
# Classification (SGMC ORIG_LABEL):
#   Carbonate:     Ked (Edwards Limestone), Kgr (Glen Rose Limestone),
#                  Kbu (Buda Limestone), Kdv (Devils River Limestone)
#   Non-carbonate: Eqc (Queen City Sand), Ecm (Cook Mountain Formation)
#
# Excluded units and rationale:
#   Kau (Austin Chalk)  — lies along the Balcones Escarpment, not on the
#                         plateau itself; topographic and climatic confound.
#   Kwa (Walnut Clay)   — classified as "undifferentiated" in SGMC; a
#                         calcareous clay marl that is hydrologically
#                         intermediate between the two groups.
#   Mixed units         — too small and lithologically heterogeneous.
#   Quaternary units    — unconsolidated deposits; CWDmax reflects
#                         alluvial position, not bedrock geology.
#
# Expected direction (H3): Non-carbonate > Carbonate.
# Rationale: clastic substrates develop deeper soils and retain water
# longer than karstified carbonates, which drain rapidly (Heilman et al.
# 2012).
#
# Approach mirrors 03_h3_geology_eel_river.R for methodological consistency.
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

r      <- terra::rast(here("data", "edwards_plateau_9ref.tif"))
r_cwd  <- r[["cwd_max"]]
r_gaps <- r[["gap_filled_months"]]


# ---- Load and classify SGMC -------------------------------------------------

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

carbonate_units     <- c("Ked", "Kgr", "Kbu", "Kdv")
non_carbonate_units <- c("Eqc", "Ecm")

sgmc <- sf::st_read(
  dsn        = here("data-raw", "USGS_StateGeologicMapCompilation_ver1.1.gdb"),
  layer      = "SGMC_Geology",
  wkt_filter = sf::st_as_text(aoi_gdb),
  quiet      = TRUE
) |>
  sf::st_transform(crs = 5070) |>
  dplyr::filter(ORIG_LABEL %in% c(carbonate_units, non_carbonate_units)) |>
  dplyr::mutate(
    geo_class = dplyr::case_when(
      ORIG_LABEL %in% carbonate_units     ~ "Carbonate",
      ORIG_LABEL %in% non_carbonate_units ~ "Non-carbonate"
    ),
    geo_class = factor(geo_class,
                       levels = c("Non-carbonate", "Carbonate"))
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
    values = c("Non-carbonate" = "#1a9850",
               "Carbonate"     = "#d73027")
  ) +
  labs(
    x = NULL,
    y = expression(paste(CWD[max], " (mm)"))
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  here("fig", "h3_edwards_boxplot.pdf"),
  plot   = p_box,
  width  = 10,
  height = 10,
  units  = "cm"
)


# ---- Density plot -----------------------------------------------------------

p_dens <- ggplot(
  df_sample,
  aes(x = cwd_max, fill = geo_class, colour = geo_class)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(
    values = c("Non-carbonate" = "#1a9850",
               "Carbonate"     = "#d73027"),
    name   = NULL
  ) +
  scale_colour_manual(
    values = c("Non-carbonate" = "#1a9850",
               "Carbonate"     = "#d73027"),
    name   = NULL
  ) +
  labs(
    x = expression(paste(CWD[max], " (mm)")),
    y = "Density"
  ) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.8))

ggsave(
  here("fig", "h3_edwards_density.pdf"),
  plot   = p_dens,
  width  = 12,
  height = 8,
  units  = "cm"
)

message("\nFigures saved to fig/")