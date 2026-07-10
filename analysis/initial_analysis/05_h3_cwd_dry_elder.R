# 05_h3_cwd_dry_elder.R
# Compares CWDmax at Elder Creek (Coastal Belt) and Dry Creek (Central Belt)
# in the Eel River Critical Zone Observatory.
#
# Outputs:
#   - data/h3_cwd_validation.csv
#   - fig/h3_cwd_validation.pdf / .png


# ---- Setup ------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(here)

stopifnot(
  "CWD raster not found. Place eel_3ref.tif in data/." =
    file.exists(here("data", "eel_3ref.tif"))
)

# Wong-palette colours: same assignment as 04_h3_geology_eel_river.R so that
# the Coastal Belt / Central Belt distinction is encoded identically across
# figures. Elder Creek sits in the Coastal Belt (TK), Dry Creek in the
# Central Belt (KJf).
site_cols <- c("Elder Creek" = "#009E73",
               "Dry Creek"   = "#D55E00")

buffer_m <- 300   # ≈ 10 pixels at 30 m resolution


# ---- Load raster ------------------------------------------------------------

r_cwd <- terra::rast(here("data", "eel_3ref.tif"))[["cwd_max"]]


# ---- Define reference sites (WGS84 / EPSG:4326) ----------------------------

# Catchment mouth coordinates from Hahm et al. (2019, WRR) Table 1.
sites <- data.frame(
  site = c("Elder Creek", "Dry Creek"),
  lon  = c(-123.6477, -123.4642),
  lat  = c(  39.7284,   39.5754)
)

pts <- sf::st_as_sf(sites, coords = c("lon", "lat"), crs = 4326) |>
  sf::st_transform(crs = 5070) |>
  sf::st_buffer(dist = buffer_m)


# ---- Extract pixels ---------------------------------------------------------

df <- terra::extract(r_cwd, terra::vect(pts), fun = NULL, na.rm = TRUE) |>
  dplyr::mutate(
    site = factor(sites$site[ID], levels = c("Elder Creek", "Dry Creek"))
  ) |>
  dplyr::filter(!is.na(cwd_max))

readr::write_csv(df, here("data", "h3_cwd_validation.csv"))


# ---- Boxplot ----------------------------------------------------------------

p <- ggplot(df, aes(x = site, y = cwd_max, fill = site)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
  scale_fill_manual(values = site_cols) +
  labs(
    x     = NULL,
    y     = expression(CWD[max]~"(mm)"),
    title = "CWD_max at Elder Creek and Dry Creek"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  here("fig", "h3_elder_dry.pdf"),
  plot   = p,
  width  = 12,
  height = 10,
  units  = "cm"
)
ggsave(
  here("fig", "h3_elder_dry.png"),
  plot   = p,
  width  = 12,
  height = 10,
  units  = "cm",
  dpi    = 600
)

message("Outputs saved to data/ and fig/")