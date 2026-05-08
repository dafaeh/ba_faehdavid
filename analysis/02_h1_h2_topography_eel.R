# 02_h1_h2_topography_eel_river.R
# Tests H1 (CWDmax decreases with elevation) and H2 (valley bottoms show
# higher CWDmax than ridges) for the Eel River focus region.
#
# H1: CWDmax decreases with elevation, driven by higher precipitation and
# lower ET at high elevations.
# H2: Valley bottoms (high TWI) exhibit higher CWDmax than ridges (low TWI),
# because lateral subsurface flow subsidies sustain higher transpiration
# deeper into the dry season (Tai et al. 2020, WRR).
#
# Inputs (data/):
#   eel_3ref.tif              CWDmax from GEE (bands: cwd_max,
#                             gap_filled_months), EPSG:5070, 30 m
#   elevation_eel.tif         FABDEM from GEE, EPSG:5070, 30 m
#
# Intermediates (data/):
#   twi_eel.tif                   TWI computed by whitebox (cached)
#   h1_h2_eel_pixels.csv          Stratified subsample used for all models
#
# Outputs (fig/):
#   h1_elevation_regression_eel.pdf
#   h2a_twi_regression_eel.pdf
#   h2b_valley_ridge_boxplot_eel.pdf
#   h2c_twi_quantile_boxplot_eel.pdf
#
# References:
#   Beven & Kirkby (1979) — TWI formula
#   Tarboton (1997)       — D-infinity flow accumulation algorithm
#   Lindsay (2016)        — WhiteboxTools
#   Hawker et al. (2022)  — FABDEM


# ---- Setup ------------------------------------------------------------------

library(terra)
library(whitebox)
library(dplyr)
library(ggplot2)
library(readr)
library(broom)
library(here)

set.seed(42)

region      <- "eel"
n_sample    <- 100000   # total pixels; divided equally across elevation bands
n_elev_band <- 10       # number of elevation bands for stratified sampling

path_cwd    <- here("data", "eel_3ref.tif")
path_elev   <- here("data", "elevation_eel.tif")
path_twi    <- here("data", "twi_eel.tif")
path_pixels <- here("data", "h1_h2_eel_pixels.csv")

dir.create(here("data", "tmp"), showWarnings = FALSE)
dir.create(here("fig"),         showWarnings = FALSE)


# ---- TWI computation (cached) -----------------------------------------------

# TWI is computed once and saved to disk. On re-runs the cached file is
# loaded instead. Workflow follows Lindsay (2016):
#   breach depressions → fill residuals → D-inf SCA → slope → TWI
# D-infinity (Tarboton 1997) is preferred over D8 because it distributes
# flow across multiple neighbours and produces smoother SCA estimates.

if (file.exists(path_twi)) {
  message("TWI cache found — skipping computation: ", path_twi)
} else {
  message("Computing TWI for ", region, " ...")
  whitebox::wbt_init()
  
  path_breached <- here("data", "tmp",
                        paste0("elevation_", region, "_breached.tif"))
  path_filled   <- here("data", "tmp",
                        paste0("elevation_", region, "_conditioned.tif"))
  path_sca      <- here("data", "tmp", paste0("sca_",   region, ".tif"))
  path_slope    <- here("data", "tmp", paste0("slope_", region, ".tif"))
  
  whitebox::wbt_breach_depressions_least_cost(
    dem    = path_elev,
    output = path_breached,
    dist   = 5,
    fill   = TRUE
  )
  
  whitebox::wbt_fill_depressions_wang_and_liu(
    dem    = path_breached,
    output = path_filled
  )
  
  whitebox::wbt_d_inf_flow_accumulation(
    input    = path_filled,
    output   = path_sca,
    out_type = "Specific Contributing Area",
    log      = FALSE
  )
  
  whitebox::wbt_slope(
    dem    = path_filled,
    output = path_slope,
    units  = "degrees"
  )
  
  whitebox::wbt_wetness_index(
    sca    = path_sca,
    slope  = path_slope,
    output = path_twi
  )
  
  message("TWI saved to: ", path_twi)
}


# ---- Stack & stratified subsample -------------------------------------------

# Elevation bands are unequally distributed across the landscape; a simple
# random sample would undersample high-elevation pixels. Equal-sized strata
# based on elevation quantiles ensures all elevation ranges are represented.
# The stratification is only for sampling — models use the full continuous
# variables, not the strata.

message("Loading rasters and building stack ...")

cwd  <- terra::rast(path_cwd)[["cwd_max"]]
elev <- terra::rast(path_elev)
twi  <- terra::rast(path_twi)

# Both elevation and TWI are already in EPSG:5070 (GEE export); resample
# to the CWDmax grid in case minor pixel-alignment differences exist.
elev <- terra::resample(elev, cwd, method = "bilinear")
twi  <- terra::resample(twi,  cwd, method = "bilinear")

stack        <- c(cwd, elev, twi)
names(stack) <- c("cwd_max", "elevation", "twi")

message("Building stratified sample (", n_elev_band, " elevation bands, ",
        n_sample / n_elev_band, " pixels each) ...")

# method = "regular" avoids random disk-seek overhead that stalls for large
# regions with high stack-NA fractions. Oversampling ensures enough valid
# pixels per elevation band after NA removal.
df_full <- terra::spatSample(
  stack,
  size   = n_sample * 3,
  method = "regular",
  na.rm  = TRUE,
  as.df  = TRUE
)

# Assign elevation bands and draw equal-sized samples per band.
n_per_band <- n_sample / n_elev_band

df <- df_full |>
  dplyr::mutate(
    elev_band = cut(
      elevation,
      breaks         = quantile(elevation, probs = seq(0, 1, 1 / n_elev_band)),
      labels         = FALSE,
      include.lowest = TRUE
    )
  ) |>
  dplyr::group_by(elev_band) |>
  dplyr::slice_sample(n = n_per_band) |>
  dplyr::ungroup() |>
  dplyr::select(-elev_band)

message(sprintf("Final sample size: %d pixels", nrow(df)))

readr::write_csv(df, path_pixels)
message("Pixel sample saved to: ", path_pixels)


# ---- H1: CWDmax ~ Elevation (continuous) ------------------------------------

message("\n=== H1: CWDmax ~ Elevation ===")

mod_elev <- lm(cwd_max ~ elevation, data = df)

print(broom::glance(mod_elev))
print(broom::tidy(mod_elev, conf.int = TRUE))

p_h1 <- ggplot(
  df,
  aes(x = elevation, y = cwd_max)) +
  geom_hex(bins = 60) +
  geom_smooth(method = "lm", colour = "firebrick", linewidth = 0.8) +
  scale_fill_viridis_c(name = "Count") +
  labs(
    x       = "Elevation (m a.s.l.)",
    y       = expression(CWD[max]~"(mm)"),
    caption = sprintf(
      "lm: slope = %.3f mm/m, R\u00b2 = %.3f, n = %d",
      broom::tidy(mod_elev)$estimate[2],
      broom::glance(mod_elev)$r.squared,
      nrow(df)
    )
  ) +
  theme_classic()

ggsave(
  filename = here("fig", "h1_elevation_regression_eel.pdf"),
  plot     = p_h1,
  width    = 16,
  height   = 12,
  units    = "cm"
)


# ---- H2a: CWDmax ~ TWI (continuous) -----------------------------------------

message("\n=== H2a: CWDmax ~ TWI ===")

mod_twi <- lm(cwd_max ~ twi, data = df)

print(broom::glance(mod_twi))
print(broom::tidy(mod_twi, conf.int = TRUE))

p_h2a <- ggplot(
  df,
  aes(x = twi, y = cwd_max)) +
  geom_hex(bins = 60) +
  geom_smooth(method = "lm", colour = "firebrick", linewidth = 0.8) +
  scale_fill_viridis_c(name = "Count") +
  labs(
    x       = "TWI",
    y       = expression(CWD[max]~"(mm)"),
    caption = sprintf(
      "lm: slope = %.3f mm/TWI unit, R\u00b2 = %.3f, n = %d",
      broom::tidy(mod_twi)$estimate[2],
      broom::glance(mod_twi)$r.squared,
      nrow(df)
    )
  ) +
  theme_classic()

ggsave(
  filename = here("fig", "h2a_twi_regression_eel.pdf"),
  plot     = p_h2a,
  width    = 16,
  height   = 12,
  units    = "cm"
)


# ---- H2b: Valley vs. Ridge (categorical) ------------------------------------

message("\n=== H2b: Valley vs. Ridge ===")

q25 <- quantile(df$twi, 0.25)
q75 <- quantile(df$twi, 0.75)

# Pixels in the middle 50% of the TWI distribution are excluded to sharpen
# the contrast between convergent valley positions and divergent ridges.
df_h2b <- df |>
  dplyr::filter(twi < q25 | twi > q75) |>
  dplyr::mutate(
    topo_class = dplyr::if_else(twi > q75, "Valley", "Ridge"),
    topo_class = factor(topo_class, levels = c("Valley", "Ridge"))
  )

message(sprintf(
  "Valley (TWI > %.2f): n = %d | Ridge (TWI < %.2f): n = %d",
  q75, sum(df_h2b$topo_class == "Valley"),
  q25, sum(df_h2b$topo_class == "Ridge")
))

p_h2b <- ggplot(
  df_h2b,
  aes(x = topo_class, y = cwd_max)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.2) +
  labs(
    x = "Topographic class",
    y = expression(CWD[max]~"(mm)")
  ) +
  theme_classic()

ggsave(
  filename = here("fig", "h2b_valley_ridge_boxplot_eel.pdf"),
  plot     = p_h2b,
  width    = 12,
  height   = 12,
  units    = "cm"
)


# ---- H2c: TWI quantile classes (visual gradient) ----------------------------

message("\n=== H2c: CWDmax by TWI quartile ===")

df_h2c <- df |>
  dplyr::mutate(
    twi_q = cut(
      twi,
      breaks         = quantile(twi, probs = 0:4 / 4),
      labels         = c("Q1 (dry)", "Q2", "Q3", "Q4 (wet)"),
      include.lowest = TRUE
    )
  )

df_h2c |>
  dplyr::group_by(twi_q) |>
  dplyr::summarise(
    n        = dplyr::n(),
    mean_mm  = mean(cwd_max),
    sd_mm    = sd(cwd_max),
    .groups  = "drop"
  ) |>
  print()

p_h2c <- ggplot(
  df_h2c,
  aes(x = twi_q, y = cwd_max)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.2) +
  labs(
    x = "TWI quartile",
    y = expression(CWD[max]~"(mm)")
  ) +
  theme_classic()

ggsave(
  filename = here("fig", "h2c_twi_quantile_boxplot_eel.pdf"),
  plot     = p_h2c,
  width    = 14,
  height   = 12,
  units    = "cm"
)

message("\nAll outputs saved to fig/ and data/.")