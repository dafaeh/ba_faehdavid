# 04_h1_h2_topography_appalachia.R
# This Script tests H1 (CWDmax decreases with elevation) and H2 (valley bottoms show
# higher CWDmax than ridges) for the Appalachian focus region. It works in the
# exact same manner as the 03_h1_h2_topography_northern_california.R script.

# ---- Setup ------------------------------------------------------------------
library(terra)
library(dplyr)
library(ggplot2)
library(readr)
library(broom)
library(here)
library(mgcv)

source(here("R", "annotate_lm.R"))
source(here("R", "lanid_filter.R"))
source(here("R", "save_fig.R"))

set.seed(42) 

# Load data
path_cwd    <- here("data", "cwd_appalachia.tif")
path_stack  <- here("data", "stack_appalachia.tif")

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

# ---- Subsample ---------------------------------------------------------------
df_full <- terra::spatSample(
  stack,
  size   = n_sample * 2,
  method = "regular",
  na.rm  = TRUE,
  xy     = TRUE,
  as.df  = TRUE
)

# Exclude irrigated pixels and pixels without LANID coverage
df_full$lanid <- terra::extract(lanid, df_full[, c("x", "y")])$lanid
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

# Add R²/n/slope to plot
plot_elev_hex <- annotate_lm(plot_elev_hex, mod_h1, slope_unit = "mm/m")

# Save plot
save_fig(plot_elev_hex, "h1_elev_hex_appalachia")

# ---- Test H2: CWDmax vs. TWI ------------------------------------------------
mod_h2 <- lm(cwd_max ~ twi, data = df)
print(broom::tidy(mod_h2, conf.int = TRUE))
print(broom::glance(mod_h2)[c("r.squared", "adj.r.squared")])

# Hexbin plot of CWDmax vs. TWI with linear and GAM trend line
plot_cwd_twi <- ggplot(df, aes(x = twi, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(name = "n Pixel", trans = "log10") +
  geom_smooth(method = "lm", colour = "firebrick", linewidth = 0.9, se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x), colour = "black",
              linewidth = 0.9, linetype = "dashed") +
  labs(
    x     = "TWI",
    y       = expression(CWD[max]~"[mm]")) +
  theme_classic()

# Add R²/n/slope. TWI is dimensionless, so the slope is reported
# per TWI unit rather than a physical unit.
plot_cwd_twi <- annotate_lm(plot_cwd_twi, mod_h2, slope_unit = "mm/TWI unit")

# Save plot
save_fig(plot_cwd_twi, "h2_twi_hex_appalachia")
