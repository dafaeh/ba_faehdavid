# 03_h1_h2_topography_eel.R
# This Script tests H1 (CWDmax decreases with elevation) and H2 (valley bottoms show
# higher CWDmax than ridges) for the Eel River focus region.
#
# H1 (elevation) uses the full Eel River focus region.
# H2 (TWI) is restricted to eel_shading_subset, used to test sensitivity to
# shading effects on DisALEXI ET (see Risks and Contingency in the proposal).

# ---- Setup ------------------------------------------------------------------
library(terra)
library(dplyr)
library(ggplot2)
library(broom)
library(here)
library(mgcv)

# annotate_lm() adds an R²/n/slope label to scatter/hexbin plots with a
# trend line (see R/annotate_lm.R for the definition).
# drop_irrigated() applies the project-wide LANID rule (see R/lanid_filter.R).
# save_fig() writes each figure to fig/ as PDF and PNG (see R/save_fig.R).
source(here("R", "annotate_lm.R"))
source(here("R", "lanid_filter.R"))
source(here("R", "save_fig.R"))

set.seed(42) # Ensure reproducibility

# Load data
path_cwd    <- here("data", "eel_9ref.tif")
path_stack  <- here("data", "stack_eel.tif")

# Define sample size
n_sample    <- 1000000

# Load CWD raster
cwd <- terra::rast(path_cwd)[["cwd_max"]]

# Load elevation and TWI from the preprocessed stack
stack_pre <- terra::rast(path_stack)
elev  <- stack_pre[["elevation"]]
twi   <- stack_pre[["twi"]]
lanid <- stack_pre[["lanid"]]

stack        <- c(cwd, elev, twi)
names(stack) <- c("cwd_max", "elevation", "twi")

# ---- Subsample (full Eel River region, used for H1) -------------------------
df_full <- terra::spatSample(
  stack,
  size   = n_sample * 2,
  method = "regular",
  na.rm  = TRUE,
  xy     = TRUE,
  as.df  = TRUE
)

# Exclude irrigated pixels and pixels without LANID coverage
df_full$lanid <- terra::extract(lanid, as.matrix(df_full[, c("x", "y")]))$lanid
df_full       <- drop_irrigated(df_full)

# min() guards against the LANID filter leaving fewer than n_sample rows,
# which slice_sample() would treat as an error.
df <- df_full |>
  dplyr::slice_sample(n = min(n_sample, nrow(df_full)))

# ---- Test H1: CWDmax vs. Elevation ------------------------------------------
mod_h1 <- lm(cwd_max ~ elevation, data = df)
print(broom::tidy(mod_h1, conf.int = TRUE))
print(broom::glance(mod_h1)[c("r.squared", "adj.r.squared")])

# Hexbin plot of CWDmax vs. elevation with linear and GAM trend line
plot_elev_hex <- ggplot(df, aes(x = elevation, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(name = "n Pixel", trans = "log10") +
  geom_smooth(method = "lm", colour = "firebrick", linewidth = 0.9) +
  geom_smooth(method = "gam", formula = y ~ s(x), colour = "black",
              linewidth = 0.9, linetype = "dashed") +
  labs(
    x       = "Elevation [m]",
    y       = expression(CWD[max]~"[mm]"))+
  theme_classic()

# Add R²/n/slope label. Slope unit is mm CWDmax per m elevation.
plot_elev_hex <- annotate_lm(plot_elev_hex, mod_h1, slope_unit = "mm/m")

# Save plot
save_fig(plot_elev_hex, "h1_elev_hex_eel")

# ---- Restrict to eel_shading_subset (used for H2) ----------------------------
# Subset defined directly in EPSG:5070 (the CRS of the rasters), NOT as a
# WGS84 rectangle. A WGS84 rectangle is a tilted quadrilateral in Albers,
# and terra::crop() crops to a vector's *bounding box*, never to the
# polygon itself -- so cropping with the transformed WGS84 rectangle would
# take the axis-aligned box around the tilted shape, roughly twice the
# intended area (the subset is narrow and tall, which maximises the effect).
# These four numbers must stay identical to ext_subset in
# 05_shading_artefact_eel.R and to subset_bbox in plot_focus_regions.R.
ext_subset <- terra::ext(
  -2336730, -2229380,   # xmin, xmax
  2138000,  2366160    # ymin, ymax
)

cwd_subset       <- terra::crop(cwd, ext_subset)
stack_pre_subset <- terra::crop(stack_pre, ext_subset)
lanid_subset     <- stack_pre_subset[["lanid"]]

stack_subset        <- c(
  cwd_subset,
  stack_pre_subset[["elevation"]],
  stack_pre_subset[["twi"]]
)
names(stack_subset) <- c("cwd_max", "elevation", "twi")

# ---- Subsample (eel_shading_subset, used for H2) -----------------------------
df_subset_full <- terra::spatSample(
  stack_subset,
  size   = n_sample * 2,
  method = "regular",
  na.rm  = TRUE,
  xy     = TRUE,
  as.df  = TRUE
)

# Exclude irrigated pixels and pixels without LANID coverage
df_subset_full$lanid <- terra::extract(
  lanid_subset, as.matrix(df_subset_full[, c("x", "y")])
)$lanid
df_subset_full       <- drop_irrigated(df_subset_full)

df_subset <- df_subset_full |>
  dplyr::slice_sample(n = min(n_sample, nrow(df_subset_full)))

# ---- Test H2: CWDmax vs. TWI ------------------------------------------------
mod_h2 <- lm(cwd_max ~ twi, data = df_subset)
print(broom::tidy(mod_h2, conf.int = TRUE))
print(broom::glance(mod_h2)[c("r.squared", "adj.r.squared")])

# Hexbin plot of CWDmax vs. TWI with linear and GAM trend line
plot_cwd_twi <- ggplot(df_subset, aes(x = twi, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(name = "n Pixel", trans = "log10") +
  geom_smooth(method = "lm", colour = "firebrick", linewidth = 0.9) +
  geom_smooth(method = "gam", formula = y ~ s(x), colour = "black",
              linewidth = 0.9, linetype = "dashed") +
  labs(
    x     = "TWI",
    y       = expression(CWD[max]~"[mm]")) +
  theme_classic()

# Add R²/n/slope label. TWI is dimensionless, so the slope is reported
# per TWI unit.
plot_cwd_twi <- annotate_lm(plot_cwd_twi, mod_h2, slope_unit = "mm/TWI unit")

# Save plot
save_fig(plot_cwd_twi, "h2_twi_hex_eel_subset")