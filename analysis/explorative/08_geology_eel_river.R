# 08_geology_eel_river.R
# Tests whether the Coastal Belt vs. Central Belt CWDmax contrast, validated
# at the point scale in 07_h3_cwd_dry_elder.R, holds across the wider Eel
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

# add_boxplot_n() adds per-group n to the boxplot below.
# drop_irrigated() applies the project-wide LANID rule (see R/lanid_filter.R).
# save_fig() writes each figure to fig/ as PDF and PNG (see R/save_fig.R).
source(here("R", "add_boxplot_n.R"))
source(here("R", "lanid_filter.R"))
source(here("R", "save_fig.R"))

set.seed(42)

# Wong-palette colours for the Coastal Belt / Central Belt distinction.
# 07_h3_cwd_dry_elder.R currently uses a single grey fill instead, so this
# encoding is not yet shared across the two H3 figures.
geo_cols <- c("Coastal Belt (TK)"  = "#009E73",
              "Central Belt (KJf)" = "#D55E00")


# Load CWD and LANID rasters
r      <- terra::rast(here("data", "eel_9ref.tif"))
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

# Exclude irrigated pixels and pixels without LANID coverage
df_pixels <- drop_irrigated(df_pixels)


# Sample
n_sample <- 100000

# min() guards against a geo_class holding fewer than n_sample pixels after
# the LANID filter, which slice_sample() would treat as an error.
df_sample <- df_pixels |>
  dplyr::group_by(geo_class) |>
  dplyr::slice_sample(n = min(n_sample, dplyr::n())) |>
  dplyr::ungroup()


# ---- Boxplot ----------------------------------------------------------------

p_box <- ggplot(
  df_sample,
  aes(x = geo_class, y = cwd_max, fill = geo_class)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
  scale_fill_manual(values = geo_cols) +
  labs(
    x = NULL,
    y = expression(CWD[max]~"(mm)")
  ) +
  theme_classic() +
  theme(legend.position = "none")

# Add per-group n label.
p_box <- p_box + add_boxplot_n(df_sample, "geo_class", "cwd_max")

save_fig(p_box, "h3_eel_river_boxplot", width = 12, height = 10)