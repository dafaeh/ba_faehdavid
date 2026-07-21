# 04_h1_h2_topography_appalachia.R
# This Script tests H1 (CWDmax decreases with elevation) and H2 (valley bottoms show
# higher CWDmax than ridges) for the Appalachian focus region. It works in the
# exact same manner as the 03_h1_h2_topography_eel.R script.

# ---- Setup ------------------------------------------------------------------
library(terra)
library(dplyr)
library(ggplot2)
library(readr)
library(broom)
library(here)
library(mgcv)

set.seed(42) # Ensure reproducibility

# Load data
path_cwd    <- here("data", "appalachia_9ref.tif")
path_stack  <- here("data", "stack_appalachia.tif")
path_pixels <- here("data", "h1_h2_appalachia_pixels.csv")

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

# Exclude irrigated pixels
df_full$lanid <- terra::extract(lanid, df_full[, c("x", "y")])$lanid
n_irrigated   <- sum(df_full$lanid == 1, na.rm = TRUE)
df_full       <- dplyr::filter(df_full, is.na(lanid) | lanid != 1)

df <- df_full |>
  dplyr::slice_sample(n = n_sample)

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

# Save plot
ggsave(
  filename = here("fig", "h1_elev_hex_appalachia.pdf"),
  plot     = plot_elev_hex,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h1_elev_hex_appalachia.png"),
  plot     = plot_elev_hex,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)

# ---- Test H2: CWDmax vs. TWI ------------------------------------------------
mod_h2 <- lm(cwd_max ~ twi, data = df)
print(broom::tidy(mod_h2, conf.int = TRUE))
print(broom::glance(mod_h2)[c("r.squared", "adj.r.squared")])

# Hexbin plot of CWDmax vs. TWI with linear and GAM trend line
plot_cwd_twi <- ggplot(df, aes(x = twi, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(name = "n Pixel", trans = "log10") +
  geom_smooth(method = "lm", colour = "firebrick", linewidth = 0.9) +
  geom_smooth(method = "gam", formula = y ~ s(x), colour = "black",
              linewidth = 0.9, linetype = "dashed") +
  labs(
    x     = "TWI",
    y       = expression(CWD[max]~"[mm]")) +
  theme_classic()

# Save plot
ggsave(
  filename = here("fig", "h2_twi_hex_appalachia.pdf"),
  plot     = plot_cwd_twi,
  width    = 16,
  height   = 12,
  units    = "cm"
)
ggsave(
  filename = here("fig", "h2_twi_hex_appalachia.png"),
  plot     = plot_cwd_twi,
  width    = 16,
  height   = 12,
  units    = "cm",
  dpi      = 600
)