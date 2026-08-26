library(here)
library(terra)
library(ggplot2)
library(tidyterra)

# ---- Einlesen & Mergen auf Disk ---------------------------------------------
files <- list.files(
  here("data", "conus_grid"),
  pattern    = "\\.tif$",
  full.names = TRUE
)

tiles <- lapply(files, terra::rast)
col   <- terra::sprc(tiles)

cwd <- terra::merge(
  col,
  filename  = here("data", "cwd_conus.tif"),
  overwrite = TRUE,
  datatype  = "FLT4S",
  gdal      = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=YES")
)

# Sanity-Check: Resolution, CRS (EPSG:5070), Bandnamen
print(cwd)
names(cwd)

# ---- Plot -------------------------------------------------------------------
ggplot() +
  tidyterra::geom_spatraster(data = cwd["cwd_max"]) +
  scale_fill_viridis_c(
    name     = expression(CWD[max] ~ "(mm)"),
    na.value = NA
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic()