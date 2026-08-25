# 08_h3_frio_canyon_transect.R
# H3 test in a tight paired-transect design at Frio Canyon (Garner State
# Park area, Edwards Plateau). Compares CWDmax across three geological
# units that sit within ~20 km of each other at near-identical climate
# and vegetation conditions, in the spirit of the Elder Creek / Dry Creek
# comparison at the Eel River (04, 05).
#
# Rationale: the regional Edwards Plateau tests (06 binary, 07 ANCOVA) are
# confounded by climate, terrain, and land-use gradients across the
# Plateau. Frio Canyon concentrates Devils River Limestone (plateau-cap,
# massive reef limestone, partly karstified), Glen Rose Limestone
# (canyon-flank, thinly bedded limestone with marl interbeds, higher
# matrix porosity), and Quaternary deposits (valley-bottom alluvium)
# inside a ~14 x 21 km box with effectively uniform climate and dominant
# rangeland / oak-juniper cover.
#
# SGMC units in the AOI (from diagnostic query):
#   Kdv  "Devils River Limestone"           Sedimentary, carbonate
#   Kgr  "Glen Rose Limestone"              Sedimentary, carbonate
#   Qu   "Quaternary deposit, undivided"    Unconsolidated, undifferentiated
#
# Note: both carbonates differ hydraulically. Devils River is a massive
# Cretaceous reef limestone prone to karst dissolution; Glen Rose is a
# cyclic alternation of limestone and marl with higher matrix porosity
# and less karst development. The comparison therefore tests storage
# differences WITHIN the carbonate class — finer than the binary
# carbonate/clastic split from 06.
#
# Hypotheses:
#   H3a: Glen Rose > Devils River, because the interbedded marl layers
#        create a deeper, more porous weathering profile that sustains
#        higher seasonal transpiration.
#   H3a': Quaternary alluvium >= Glen Rose, since unconsolidated valley
#        fill provides the greatest accessible rooting-zone storage.
#   Alternative: Devils River >= Glen Rose, if karst-mediated deep root
#        access through dissolution features dominates over matrix storage.
#
# Inputs:
#   data/edwards_plateau_9ref.tif        cwd_max + gap_filled_months
#   data-raw/USGS_StateGeologicMapCompilation_ver1.1.gdb   SGMC v1.1
#
# Outputs:
#   data/h3_frio_canyon_samples.csv
#   data/h3_frio_canyon_unit_inventory.csv
#   data/h3_frio_canyon_desc.csv
#   fig/h3_frio_canyon_map.pdf
#   fig/h3_frio_canyon_boxplot.pdf
#   fig/h3_frio_canyon_density.pdf


# ---- Setup ------------------------------------------------------------------

library(terra)
library(tidyterra)
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(here)

set.seed(42)

gap_max <- 6


# ---- Define AOI -------------------------------------------------------------

aoi_bbox <- c(
  xmin = -99.79, ymin = 29.55,
  xmax = -99.66, ymax = 29.74
)

aoi_wgs84 <- sf::st_as_sfc(sf::st_bbox(aoi_bbox, crs = sf::st_crs(4326)))
aoi_5070  <- sf::st_transform(aoi_wgs84, crs = 5070)


# ---- Load and crop CWD raster -----------------------------------------------

stopifnot(
  "CWD raster not found. Place edwards_plateau_9ref.tif in data/." =
    file.exists(here("data", "edwards_plateau_9ref.tif"))
)

r_full <- terra::rast(here("data", "edwards_plateau_9ref.tif"))
r      <- terra::crop(r_full, terra::vect(aoi_5070))
r_cwd  <- r[["cwd_max"]]


# ---- Load and classify SGMC -------------------------------------------------

stopifnot(
  "SGMC not found. Run 01_download_data.R first." =
    file.exists(here("data-raw", "USGS_StateGeologicMapCompilation_ver1.1.gdb"))
)

gdb_crs <- sf::st_crs(
  sf::st_read(
    here("data-raw", "USGS_StateGeologicMapCompilation_ver1.1.gdb"),
    query = "SELECT * FROM SGMC_Geology LIMIT 1",
    quiet = TRUE
  )
)
aoi_gdb <- sf::st_transform(aoi_wgs84, crs = gdb_crs)

sgmc_all <- sf::st_read(
  dsn        = here("data-raw", "USGS_StateGeologicMapCompilation_ver1.1.gdb"),
  layer      = "SGMC_Geology",
  wkt_filter = sf::st_as_text(aoi_gdb),
  quiet      = TRUE
) |>
  sf::st_transform(crs = 5070) |>
  sf::st_intersection(aoi_5070)

# Compute area on the sf object directly — avoids the terra::geom() /
# sf geometry column name conflict inside dplyr::mutate().
sgmc_all$area_km2 <- as.numeric(sf::st_area(sgmc_all)) / 1e6

inventory <- sgmc_all |>
  sf::st_drop_geometry() |>
  dplyr::group_by(ORIG_LABEL, UNIT_NAME, GENERALIZED_LITH) |>
  dplyr::summarise(
    area_km2 = sum(area_km2, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(area_km2))

readr::write_csv(
  inventory,
  here("data", "h3_frio_canyon_unit_inventory.csv")
)

message("\n=== SGMC units within the Frio Canyon AOI ===")
print(inventory, width = Inf)


# ---- Classify the three target units ----------------------------------------

# Match on ORIG_LABEL directly — the diagnostic query confirmed the exact
# codes present in this AOI.
sgmc <- sgmc_all |>
  dplyr::mutate(
    lithology = dplyr::case_when(
      ORIG_LABEL == "Kdv" ~ "Devils River Lst. (plateau cap)",
      ORIG_LABEL == "Kgr" ~ "Glen Rose Lst. (canyon flank)",
      ORIG_LABEL == "Qu"  ~ "Quaternary (valley bottom)",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(lithology)) |>
  dplyr::mutate(
    lithology = factor(
      lithology,
      levels = c(
        "Devils River Lst. (plateau cap)",
        "Glen Rose Lst. (canyon flank)",
        "Quaternary (valley bottom)"
      )
    )
  )

message("\nPolygons per target unit:")
print(dplyr::count(sf::st_drop_geometry(sgmc), lithology))


# ---- Rasterize and extract pixels -------------------------------------------

sgmc_vect <- terra::vect(sgmc)
sgmc_vect$lithology_int <- as.integer(sgmc$lithology)

litho_rast <- terra::rasterize(sgmc_vect, r_cwd, field = "lithology_int")
names(litho_rast) <- "lithology"

stacked <- c(r, litho_rast)

# Full-pixel extraction: the AOI is small enough (~250 km^2) that 30 m
# pixels stay in the hundreds of thousands.
df_all <- terra::as.data.frame(stacked, xy = TRUE, na.rm = TRUE) |>
  tibble::as_tibble() |>
  dplyr::mutate(
    lithology = factor(
      lithology,
      levels = seq_along(levels(sgmc$lithology)),
      labels = levels(sgmc$lithology)
    )
  )


# ---- Quality filter ---------------------------------------------------------

df <- df_all |>
  dplyr::filter(
    !is.na(cwd_max),
    gap_filled_months <= gap_max
  )

message("\nValid pixels per unit after filtering:")
print(dplyr::count(df, lithology))

readr::write_csv(df, here("data", "h3_frio_canyon_samples.csv"))


# ---- Descriptive statistics -------------------------------------------------

df_desc <- df |>
  dplyr::group_by(lithology) |>
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

readr::write_csv(df_desc, here("data", "h3_frio_canyon_desc.csv"))


# ---- Colour palette ---------------------------------------------------------

litho_cols <- c(
  "Devils River Lst. (plateau cap)" = "#d73027",
  "Glen Rose Lst. (canyon flank)"    = "#fdae61",
  "Quaternary (valley bottom)"       = "#4575b4"
)


# ---- Map plot ---------------------------------------------------------------

sgmc_outline <- sgmc |>
  dplyr::group_by(lithology) |>
  dplyr::summarise(.groups = "drop")

p_map <- ggplot() +
  tidyterra::geom_spatraster(data = r_cwd) +
  scale_fill_viridis_c(
    name     = expression(CWD[max]~"(mm)"),
    na.value = NA,
    option   = "G"
  ) +
  geom_sf(
    data      = sgmc_outline,
    aes(colour = lithology),
    fill      = NA,
    linewidth = 0.5
  ) +
  scale_colour_manual(values = litho_cols, name = NULL) +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    legend.position  = "right",
    legend.box       = "vertical",
    legend.spacing.y = unit(0.2, "cm")
  )

ggsave(
  here("fig", "h3_frio_canyon_map.pdf"),
  plot   = p_map,
  width  = 18,
  height = 14,
  units  = "cm"
)


# ---- Boxplot ----------------------------------------------------------------

p_box <- ggplot(df, aes(x = lithology, y = cwd_max, fill = lithology)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
  scale_fill_manual(values = litho_cols) +
  labs(
    x = NULL,
    y = expression(CWD[max]~"(mm)")
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 15, hjust = 1)
  )

ggsave(
  here("fig", "h3_frio_canyon_boxplot.pdf"),
  plot   = p_box,
  width  = 12,
  height = 10,
  units  = "cm"
)


# ---- Density plot -----------------------------------------------------------

p_dens <- ggplot(
  df,
  aes(x = cwd_max, fill = lithology, colour = lithology)
) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = litho_cols, name = NULL) +
  scale_colour_manual(values = litho_cols, name = NULL) +
  labs(
    x = expression(CWD[max]~"(mm)"),
    y = "Density"
  ) +
  theme_classic() +
  theme(legend.position = c(0.75, 0.8))

ggsave(
  here("fig", "h3_frio_canyon_density.pdf"),
  plot   = p_dens,
  width  = 14,
  height = 9,
  units  = "cm"
)


# ---- H3b: Residualvarianz nach TWI-Adjustierung -----------------------------

# Für jede Einheit separat lm(cwd_max ~ twi) fitten und den RMSE der
# Residuen berechnen. Was nach Herausrechnen der topographischen
# Wasserumverteilung übrig bleibt ist die unexplained spatial variance.
# H3b: Devils River hat höheren RMSE als Glen Rose bei vergleichbarem
# oder niedrigerem Median — Karststrukturen schaffen diskrete Hochpunkte,
# die TWI nicht erklären kann.
#
# Voraussetzung: TWI-Raster muss für das Frio-Canyon-AOI verfügbar sein.

stopifnot(
  "TWI raster not found — H3b requires TWI." =
    file.exists(here("data", "CONUS_TWI_epsg5072_30m_unmasked.tif"))
)

twi_conus <- terra::rast(here("data", "CONUS_TWI_epsg5072_30m_unmasked.tif"))

aoi_in_twi_crs <- terra::project(
  terra::vect(terra::ext(r_cwd), crs = terra::crs(r_cwd)),
  terra::crs(twi_conus)
)
twi_crop    <- terra::crop(twi_conus, aoi_in_twi_crs)
twi_aligned <- terra::project(twi_crop, r_cwd, method = "bilinear")
names(twi_aligned) <- "twi"

df_twi <- df |>
  dplyr::mutate(twi = terra::extract(twi_aligned, cbind(x, y))[, 1]) |>
  dplyr::filter(!is.na(twi))

# RMSE und Residual-SD pro Einheit.
h3b_var <- df_twi |>
  dplyr::group_by(lithology) |>
  dplyr::summarise(
    n            = dplyr::n(),
    median_mm    = median(cwd_max),
    sd_raw_mm    = sd(cwd_max),
    rmse_resid   = {
      mod_twi <- lm(cwd_max ~ twi, data = dplyr::pick(everything()))
      sqrt(mean(resid(mod_twi)^2))
    },
    r2_twi       = {
      mod_twi <- lm(cwd_max ~ twi, data = dplyr::pick(everything()))
      summary(mod_twi)$r.squared
    },
    .groups = "drop"
  )

message("\n=== H3b: Residualvarianz nach TWI-Adjustierung ===")
print(as.data.frame(h3b_var), digits = 1)

readr::write_csv(h3b_var, here("data", "h3_frio_canyon_h3b_variance.csv"))

# Plot: rohe SD vs. Residual-RMSE nebeneinander, um zu zeigen wie viel
# die TWI-Adjustierung ausmacht.
h3b_long <- h3b_var |>
  dplyr::select(lithology, sd_raw_mm, rmse_resid) |>
  tidyr::pivot_longer(
    cols      = c(sd_raw_mm, rmse_resid),
    names_to  = "measure",
    values_to = "value"
  ) |>
  dplyr::mutate(
    measure = dplyr::recode(
      measure,
      sd_raw_mm  = "Raw SD",
      rmse_resid = "Residual RMSE\n(after TWI)"
    )
  )

p_var <- ggplot(
  h3b_long,
  aes(x = lithology, y = value, fill = lithology)
) +
  geom_col() +
  scale_fill_manual(values = litho_cols) +
  facet_wrap(~ measure) +
  labs(
    x = NULL,
    y = "CWDmax variability (mm)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 15, hjust = 1)
  )

ggsave(
  here("fig", "h3_frio_canyon_h3b_variance.pdf"),
  plot   = p_var,
  width  = 16,
  height = 9,
  units  = "cm"
)

message("\nAll outputs written to data/ and fig/.")