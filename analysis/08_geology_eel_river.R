# 08_geology_eel_river.R
# Tests whether the Coastal Belt vs. Central Belt CWDmax contrast, validated
# at the point scale in 08_h3_cwd_dry_elder.R, holds across the wider Eel
# River region. Extracts CWDmax over all TK/KJf polygons in the focus region
# and tests H3 (geological substrate) at the formation scale, beyond the
# two specific catchments studied by Hahm et al. (2019).

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


# Load CWD and LANID rasters
r      <- terra::rast(here("data", "eel_3ref.tif"))
r_cwd  <- r[["cwd_max"]]
r_gaps <- r[["gap_filled_months"]]
lanid <- terra::rast(here("data", "stack_eel.tif"))[["lanid"]]


# Load and classify SGMC 
sgmc <- sf::st_read(
  here("data", "geology_eel.gpkg"),
  quiet = TRUE
) |>
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
# downstream if needed. cells = TRUE returns the grid cell number of each
# extracted pixel, used to look up LANID on the shared grid.
extracted_cwd  <- terra::extract(r_cwd,  terra::vect(sgmc), fun = NULL,
                                 cells = TRUE)
extracted_gaps <- terra::extract(r_gaps, terra::vect(sgmc), fun = NULL)

# LANID value at each extracted pixel, via cell number (same grid as r_cwd).
lanid_at <- lanid[extracted_cwd$cell][["lanid"]]

df_pixels <- dplyr::tibble(
  ID                = extracted_cwd$ID,
  cwd_max           = extracted_cwd$cwd_max,
  gap_filled_months = extracted_gaps$gap_filled_months,
  geo_class         = sgmc$geo_class[extracted_cwd$ID],
  orig_label        = sgmc$ORIG_LABEL[extracted_cwd$ID],
  lanid             = lanid_at
) |>
  dplyr::filter(!is.na(cwd_max))

# Exclude irrigated pixels
n_irrigated <- sum(df_pixels$lanid == 1, na.rm = TRUE)
df_pixels   <- dplyr::filter(df_pixels, is.na(lanid) | lanid != 1)


# Sample
n_sample <- 100000

df_sample <- df_pixels |>
  dplyr::group_by(geo_class) |>
  dplyr::slice_sample(n = n_sample) |>
  dplyr::ungroup()


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