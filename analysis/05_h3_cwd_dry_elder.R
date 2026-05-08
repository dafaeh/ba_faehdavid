# 05_cwd_smax_validation.R
# Compares CWDmax values from the GEE export (eel_3ref.tif) with
# independently measured rooting-zone water storage capacity (Smax) for
# Elder Creek (Coastal Belt) and Dry Creek (Central Belt) in the Eel
# River Critical Zone Observatory.
#
# Following Stocker et al. (2023, Nat. Geosci.), CWDmax should approximate
# Smax because CWD measures the cumulative water drawn from storage. If the
# dry season is severe enough to fully deplete the rooting zone in at least
# one study year, CWDmax converges on Smax. If not, CWDmax is a lower bound.
#
# Reference Smax values:
#   Hahm et al. (2019, GRL) — borehole measurements + groundwater monitoring
#     Elder Creek (Coastal Belt): Smax 300–400 mm
#     Dry Creek   (Central Belt): Smax 120–200 mm
#
#   Dralle et al. (2020, ERL) — vegetation-as-sensor inversion
#     Elder Creek: inferred Smax ~303 mm
#     Dry Creek:   inferred Smax ~184 mm
#
# Site coordinates are taken directly from Table 1 and the Site
# Description section of Hahm et al. (2019, GRL).
#
# Outputs:
#   - data/cwd_validation_sites.csv
#   - data/cwd_validation_summary.csv
#   - fig/cwd_validation_plot.pdf


# ---- Setup ------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(here)


# ---- Load raster ------------------------------------------------------------

r_cwd <- terra::rast(here("data", "eel_3ref.tif"))[["cwd_max"]]

# Verify CRS: all spatial operations assume EPSG:5070 (CONUS Albers).
crs_code <- terra::crs(r_cwd, describe = TRUE)$code
if (!grepl("5070", crs_code)) {
  warning("Unexpected CRS: ", crs_code, ". Expected EPSG:5070.")
}


# ---- Define reference sites (WGS84 / EPSG:4326) -----------------------------

# Catchment outlets: Hahm et al. (2019, GRL) Table 1.
# Intensively monitored slopes: Hahm et al. (2019, GRL) Site Description.
sites_wgs84 <- data.frame(
  site_id   = c("elder_mouth", "rivendell", "dry_mouth", "dry_ridge"),
  watershed = c("Elder Creek", "Elder Creek", "Dry Creek",   "Dry Creek"),
  belt      = c("Coastal Belt", "Coastal Belt", "Central Belt", "Central Belt"),
  lon       = c(-123.6477, -123.6451, -123.4642, -123.4733),
  lat       = c(  39.7284,   39.7290,   39.5754,   39.5678)
)

pts_wgs84 <- sf::st_as_sf(
  sites_wgs84,
  coords = c("lon", "lat"),
  crs    = 4326
)

# Reproject to match the GEE export CRS.
pts_5070 <- sf::st_transform(pts_wgs84, crs = 5070)


# ---- Buffer extraction ------------------------------------------------------

# Single-point extraction is sensitive to exact pixel position; buffer
# statistics across three radii are more robust for catchment-scale sites.
buffer_radii_m <- c(100, 250, 500)  # 100 m ≈ 3 pixels, 500 m ≈ 17 pixels

extract_buffer <- function(radius_m) {
  pts_buf <- sf::st_buffer(pts_5070, dist = radius_m)
  terra::extract(r_cwd, terra::vect(pts_buf), fun = NULL, na.rm = TRUE) |>
    dplyr::rename(cwd_max_mm = names(r_cwd)) |>
    dplyr::group_by(ID) |>
    dplyr::summarise(
      n_pixel       = dplyr::n(),
      cwd_median_mm = median(cwd_max_mm, na.rm = TRUE),
      cwd_mean_mm   = mean(cwd_max_mm,   na.rm = TRUE),
      cwd_sd_mm     = sd(cwd_max_mm,     na.rm = TRUE),
      cwd_q25_mm    = quantile(cwd_max_mm, 0.25, na.rm = TRUE),
      cwd_q75_mm    = quantile(cwd_max_mm, 0.75, na.rm = TRUE),
      .groups       = "drop"
    ) |>
    dplyr::mutate(
      radius_m  = radius_m,
      site_id   = pts_5070$site_id[ID],
      belt      = pts_5070$belt[ID],
      watershed = pts_5070$watershed[ID]
    )
}

results_buf <- purrr::map(buffer_radii_m, extract_buffer) |>
  dplyr::bind_rows()

message("\n=== Buffer extraction (median CWDmax per radius) ===")
results_buf |>
  dplyr::select(site_id, belt, radius_m, n_pixel, cwd_median_mm, cwd_sd_mm) |>
  print(n = Inf)


# ---- Reference Smax values --------------------------------------------------

# Smax intervals from direct field measurements (Hahm et al. 2019, GRL)
# and the vegetation-as-sensor inversion (Dralle et al. 2020, ERL).
smax_ref <- data.frame(
  watershed          = c("Elder Creek", "Dry Creek"),
  belt               = c("Coastal Belt", "Central Belt"),
  smax_field_low_mm  = c(300, 120),
  smax_field_high_mm = c(400, 200),
  smax_field_source  = c("Hahm2019_GRL", "Hahm2019_GRL"),
  smax_inv_mm        = c(303, 184),
  smax_inv_source    = c("Dralle2020_ERL", "Dralle2020_ERL")
)


# ---- Watershed aggregation --------------------------------------------------

# Aggregate both sites per catchment using the 250 m buffer as the
# primary result; compare against reference Smax values.
cwd_per_watershed <- results_buf |>
  dplyr::filter(radius_m == 250) |>
  dplyr::group_by(watershed, belt) |>
  dplyr::summarise(
    cwd_median_mm = median(cwd_median_mm, na.rm = TRUE),
    cwd_sd_mm     = median(cwd_sd_mm,     na.rm = TRUE),
    .groups       = "drop"
  ) |>
  dplyr::left_join(smax_ref, by = c("watershed", "belt"))

message("\n=== Watershed comparison: CWDmax vs. Smax reference ===")
print(cwd_per_watershed, width = Inf)


# ---- Validation plot --------------------------------------------------------

# Show CWDmax medians for all buffer radii alongside the field-measured
# Smax range (shaded band) and vegetation-inversion estimate (dashed line).
df_plot <- results_buf |>
  dplyr::group_by(watershed, belt, radius_m) |>
  dplyr::summarise(
    cwd_median_mm = median(cwd_median_mm, na.rm = TRUE),
    .groups       = "drop"
  ) |>
  dplyr::left_join(smax_ref, by = c("watershed", "belt")) |>
  dplyr::mutate(label = paste0(belt, "\n(", watershed, ")"))

p <- ggplot(df_plot) +
  
  # Smax range from borehole measurements (Hahm et al. 2019, GRL)
  geom_rect(
    aes(xmin = smax_field_low_mm, xmax = smax_field_high_mm,
        ymin = -Inf, ymax = Inf, fill = belt),
    alpha = 0.15
  ) +
  
  # Smax from vegetation-as-sensor inversion (Dralle et al. 2020, ERL)
  geom_vline(
    aes(xintercept = smax_inv_mm, color = belt),
    linetype  = "dashed",
    linewidth = 0.8
  ) +
  
  # CWDmax medians for three buffer radii
  geom_point(
    aes(x = cwd_median_mm, y = factor(radius_m),
        color = belt, shape = belt),
    size = 3
  ) +
  
  facet_wrap(~ label, ncol = 2) +
  
  scale_color_manual(
    values = c("Coastal Belt" = "#2166ac", "Central Belt" = "#d6604d"),
    name   = NULL
  ) +
  scale_fill_manual(
    values = c("Coastal Belt" = "#2166ac", "Central Belt" = "#d6604d"),
    name   = NULL
  ) +
  scale_shape_manual(
    values = c("Coastal Belt" = 16, "Central Belt" = 17),
    name   = NULL
  ) +
  
  labs(
    x        = expression(paste("CWD"[max], " or ", S[max], "  (mm)")),
    y        = "Buffer radius (m)",
    subtitle = paste(
      "Shaded band: Smax range from borehole measurements (Hahm et al. 2019, GRL).",
      "\nDashed line: Smax from vegetation inversion (Dralle et al. 2020, ERL).",
      "\nPoints: median CWDmax for three buffer radii around reference sites."
    ),
    caption  = "Site coordinates: Hahm et al. (2019, GRL) Table 1 & Site Description."
  ) +
  
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  filename = here("fig", "cwd_validation_plot.pdf"),
  plot     = p,
  width    = 18,
  height   = 10,
  units    = "cm"
)

message("Plot saved: fig/cwd_validation_plot.pdf")


# ---- Save results -----------------------------------------------------------

readr::write_csv(results_buf,       here("data", "cwd_validation_sites.csv"))
readr::write_csv(cwd_per_watershed, here("data", "cwd_validation_summary.csv"))

message("Tables saved to data/")