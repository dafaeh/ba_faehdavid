# 02_preprocess_data.R
# Builds one multi-band raster stack containing the important variables per focus region, plus a
# region-cropped geology layer, and writes them to data/. 
# Run after 01_download_data.R

# Writes per region:
#   data/stack_<region>.tif    Bands: elevation, northness, slope, twi, nlcd,
#                              lanid, on the region's CWD grid.
#   data/geology_<region>.gpkg SGMC polygons, cropped, attributes unclassified 
#                              because different classifications are needed for 
#                              different analyses.

# Grid: <region>_9ref.tif defines the template (EPSG:5070, 30 m). Every band is
# resampled onto it.

# If you want to add a region: export it from GEE, place its exports in
# data/ (<region>_9ref.tif etc.), add its name below, re-run. No coordinates
# are needed here -- every region's extent and CRS come from its exported
# <region>_9ref.tif grid template, so this script and the GEE exports can
# never drift out of sync.

# ---- Setup ------------------------------------------------------------------
library(terra)
library(sf)
library(here)

dir.create(here("data"), showWarnings = FALSE)

# NLCD 2019 COG, streamed via /vsicurl/ (only overlapping tiles with the focus regions are read).
nlcd_url <- "/vsicurl/https://storage.googleapis.com/feddata-r/nlcd/2019_Land_Cover_L48.tif"

# ---- Focus regions ----------------------------------------------------------
# Just the names. Each region's extent/CRS is read from its exported grid
# template <region>_9ref.tif, not from hardcoded coordinates.
focus_regions <- c("eel", "edwards_plateau", "appalachia", "high_plains")

# ---- Helpers ----------------------------------------------------------------

# Path to a region's CWD grid template.
template_path <- function(region) {
  here("data", paste0(region, "_9ref.tif"))
}

# Crop to the region first (template extent reprojected into the raster's
# CRS), then resample onto the template grid. Cropping first keeps the
# expensive reproject/resample small.
align_to_template <- function(raw, template, method) {
  poly     <- terra::as.polygons(terra::ext(template), crs = terra::crs(template))
  crop_ext <- terra::ext(terra::project(poly, terra::crs(raw)))
  raw_crop <- terra::crop(raw, crop_ext)
  terra::project(raw_crop, template, method = method)
}

# ---- Loop over focus regions and builds per-region stack and geology --------

for (region in focus_regions) {
  message("=== Processing region: ", region, " ===")
  
  # Define grid template
  path_template <- template_path(region)
  if (!file.exists(path_template)) {
    stop(
      "Missing grid template ", path_template, ". Export it from GEE ",
      "(main.js) and place it in data/ before running this script."
    )
  }
  template <- terra::rast(path_template)[[1]]   # any band defines the grid
  
  # ---- elevation, northness, slope --------------------------------------
  path_elev <- here("data", paste0("elevation_", region, ".tif"))
  if (!file.exists(path_elev)) {
    stop(
      "Missing ", path_elev, ". Export elevation for '", region,
      "' from GEE and place it in data/ before running this script."
    )
  }
  
  elevation <- terra::resample(terra::rast(path_elev), template,
                               method = "bilinear")
  names(elevation) <- "elevation"
  
  # northness = cos(aspect): +1 north, -1 south.
  northness <- cos(terra::terrain(elevation, "aspect", unit = "radians"))
  names(northness) <- "northness"
  
  slope <- terra::terrain(elevation, "slope", unit = "degrees")
  names(slope) <- "slope"
  
  # ---- twi (continuous, from data-raw) ---------------------------------
  # The Zenodo TWI file has no CRS in its header, though the data is
  # EPSG:5072. Assign it before projecting (this labels metadata only, moves
  # no pixels). Without it project() errors or misaligns TWI by the small
  # 5072-vs-5070 datum offset
  path_twi_raw <- here("data-raw", "CONUS_TWI_epsg5072_30m_unmasked.tif")
  if (!file.exists(path_twi_raw)) {
    stop("Missing ", path_twi_raw, ". Run 01_download_data.R first.")
  }
  twi_raw <- terra::rast(path_twi_raw)
  if (is.na(terra::crs(twi_raw, describe = TRUE)$code)) {
    terra::crs(twi_raw) <- "EPSG:5072"
  }
  twi <- align_to_template(twi_raw, template, method = "bilinear")
  names(twi) <- "twi"
  
  # ---- nlcd (categorical, streamed, kept raw) ---------------------------
  # Kept raw (no reclassification). levels() <- NULL strips the factor RAT,
  # otherwise names() reports the RAT column. Overwriting the band name on
  # write and breaking stack_pre[["nlcd"]] downstream. Cell values unchanged.
  raw_nlcd <- terra::rast(nlcd_url)
  levels(raw_nlcd) <- NULL
  nlcd <- align_to_template(raw_nlcd, template, method = "near")
  names(nlcd) <- "nlcd"
  
  # ---- lanid (categorical, from data-raw) ------------------------------
  path_lanid <- here("data-raw", "lanid2017.tif")
  if (!file.exists(path_lanid)) {
    stop("Missing ", path_lanid, ". Run 01_download_data.R first.")
  }
  lanid <- align_to_template(terra::rast(path_lanid), template, method = "near")
  names(lanid) <- "lanid"
  
  # ---- Assemble and write the stack -------------------------------------
  # Continuous and categorical variables share Float32. Integer codes (e.g. 82.0) 
  # are not affected by this.
  stack <- c(elevation, northness, slope, twi, nlcd, lanid)
  expected_names <- c("elevation", "northness", "slope", "twi", "nlcd", "lanid")
  
  # Catch silent band-name corruption here (see nlcd note).
  stopifnot(
    "Band names in stack do not match expected names" =
      identical(names(stack), expected_names)
  )
  
  path_stack <- here("data", paste0("stack_", region, ".tif"))
  terra::writeRaster(
    stack, path_stack,
    overwrite = TRUE,
    datatype  = "FLT4S",
    gdal      = c("COMPRESS=DEFLATE", "TILED=YES")
  )
  message("  Wrote ", path_stack, " (bands: ",
          paste(names(stack), collapse = ", "), ")")
  
  # ---- geology (vector, cropped, unclassified) --------------------------
  path_gdb <- here("data-raw", "USGS_StateGeologicMapCompilation_ver1.1.gdb")
  if (!file.exists(path_gdb)) {
    warning("Missing ", path_gdb, " -- skipping geology for '", region,
            "'. Run 01_download_data.R first.")
  } else {
    # LIMIT 1: only need the CRS, not the actual data.
    gdb_crs <- sf::st_crs(sf::st_read(
      path_gdb,
      query = "SELECT * FROM SGMC_Geology LIMIT 1",
      quiet = TRUE
    ))
    
    # Filter from the template extent (the 5070 grid)
    aoi_gdb <- sf::st_transform(
      sf::st_as_sfc(sf::st_bbox(template)),
      crs = gdb_crs
    )
    
    # wkt_filter: st_read scans only overlapping polygons, not all of CONUS.
    geology <- sf::st_read(
      dsn        = path_gdb,
      layer      = "SGMC_Geology",
      wkt_filter = sf::st_as_text(aoi_gdb),
      quiet      = TRUE
    ) |>
      sf::st_transform(crs = 5070)
    
    path_geo <- here("data", paste0("geology_", region, ".gpkg"))
    sf::st_write(geology, path_geo, delete_dsn = TRUE, quiet = TRUE)
    message("  Wrote ", path_geo, " (", nrow(geology), " polygons)")
  }
  
  message("Done: ", region)
}

message("\nAll regions processed. Stacks and geology written to data/.")