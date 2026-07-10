# 02_preprocess_data.R
# Builds one multi-band raster stack containing the important variables per focus region, plus a
# region-cropped geology layer, and writes them to data/. This makes the analysis 
# scripts downstream more compact, as they only need to read these outputs. 
# Run after 01_download_data.R

# Writes per region:
#   data/stack_<region>.tif    Bands: elevation, northness, slope, twi, nlcd,
#                              lanid, on the region's CWD grid.
#   data/geology_<region>.gpkg SGMC polygons, cropped, attributes unclassified 
#                              because different classifications are needed for 
#                              different analyses.

# Reads: GEE exports <region>_9ref.tif (CWD, used as
# the grid template) and elevation_<region>.tif. Also TWI, LANID, SGMC from
# data-raw (see script 01). NLCD is streamed from a COG.

# Grid: <region>_9ref.tif defines the template (EPSG:5070, 30 m). Every band is
# resampled onto it.

# If you want to add a region: add its bounding box to focus_regions (keep in sync with config.js
# SITES by hand), place its GEE exports in data, re-run.

# ---- Setup ------------------------------------------------------------------
library(terra)
library(sf)
library(here)

dir.create(here("data"), showWarnings = FALSE)

# NLCD 2019 COG, streamed via /vsicurl/ (only overlapping tiles with the focus regions are read).
nlcd_url <- "/vsicurl/https://storage.googleapis.com/feddata-r/nlcd/2019_Land_Cover_L48.tif"

# ---- Focus region bounding boxes --------------------------------------------
# Values and order (xmin, ymin, xmax, ymax) match the ee.Geometry.Rectangle()
# arrays in config.js. Copy-paste the four numbers for each region from there into this script.
focus_regions <- list(
  eel             = c(-123.74888168761035, 39.29196540299779,
                      -121.09019028136035, 41.29344890598067),
  
  edwards_plateau = c(-100.27776465680033, 28.789320221427214,
                      -96.21282325055033, 32.660183285641914),
  
  appalachia      = c( -84.05115497183012, 37.16650782156939,
                       -78.63489520620512, 40.97746916994725),
  
  high_plains     = c(-103.73661176242638, 36.326001518952935,
                      -97.38102094211388, 40.588750552430916)
)

# ---- Helpers ----------------------------------------------------------------

# Path to a region's CWD grid template.
template_path <- function(region) {
  here("data", paste0(region, "_9ref.tif"))
}

# Crop a raster to the region, then resample onto the template grid.
# Cropping first (via the template extent reprojected into the raster's CRS)
# This keeps the expensive reproject and resample step small.
align_to_template <- function(raw, template, method) {
  poly     <- terra::as.polygons(terra::ext(template), crs = terra::crs(template))
  crop_ext <- terra::ext(terra::project(poly, terra::crs(raw)))
  raw_crop <- terra::crop(raw, crop_ext)
  terra::project(raw_crop, template, method = method)
}

# ---- Loop over focus regions and builds per-region stack and geology --------

for (region in names(focus_regions)) {
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
  # EPSG:5072 -> aligned to the 5070 template in one step. Bilinear.
  path_twi_raw <- here("data-raw", "CONUS_TWI_epsg5072_30m_unmasked.tif")
  if (!file.exists(path_twi_raw)) {
    stop("Missing ", path_twi_raw, ". Run 01_download_data.R first.")
  }
  twi <- align_to_template(terra::rast(path_twi_raw), template,
                           method = "bilinear")
  names(twi) <- "twi"
  
    # ---- nlcd (categorical, streamed, kept raw) ---------------------------
    # Native codes for the vegetation types kept (no reclassification).
    # levels() <- NULL strips the factor RAT: otherwise names() reports the RAT
    # column ("NLCD Land Cover Class"), which overwrites the band name on write
    # and breaks stack_pre[["nlcd"]] downstream. Cell values are unchanged.
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
    
    # st_bbox() needs names. focus_regions stores plain (xmin,ymin,xmax,ymax)
    # vectors to match GEE's Rectangle() order, so name them here.
    bb <- stats::setNames(focus_regions[[region]], c("xmin", "ymin", "xmax", "ymax"))
    aoi_gdb <- sf::st_transform(
      sf::st_as_sfc(sf::st_bbox(bb, crs = sf::st_crs(4326))),
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