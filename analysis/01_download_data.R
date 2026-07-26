# 01_download_data.R
# Downloads SGMC, LANID and TWI into data-raw. Run once before anything else.
#
# Not downloaded here:
#   NLCD: streamed and baked into each stack by 02_preprocess_data.R
#   CWD (<region>_9ref.tif) and elevation (elevation_<region>.tif): GEE
#     exports, place manually in data/.

# ---- Setup ------------------------------------------------------------------
library(here)

# Create data-raw folder in case it does not exist yet
dir.create(here("data-raw"), showWarnings = FALSE)

# 60s default is too short for SGMC/LANID. TWI (~41 GB) uses curl.
options(timeout = 600)

# ---- SGMC (geology) ---------------------------------------------------------
path_sgmc_zip <- here("data-raw", "USGS_SGMC_Geodatabase.zip")
path_sgmc_gdb <- here("data-raw", "USGS_StateGeologicMapCompilation_ver1.1.gdb")

if (!file.exists(path_sgmc_gdb)) {
  message("Downloading SGMC geodatabase (~400 MB) ...")
  download.file(
    url      = paste0(
      "https://www.sciencebase.gov/catalog/file/get/",
      "5888bf4fe4b05ccb964bab9d",
      "?f=__disk__24/f6/e1/",
      "24f6e139c181fd9fe43df2aaf7f50b1c5b3b6297"
    ),
    destfile = path_sgmc_zip,
    mode     = "wb"
  )
  
  message("Unzipping ...")
  unzip(path_sgmc_zip, exdir = here("data-raw"))
  file.rename(
    here("data-raw", "USGS_SGMC_Geodatabase",
         "USGS_StateGeologicMapCompilation_ver1.1.gdb"),
    path_sgmc_gdb
  )
  
  unlink(here("data-raw", "USGS_SGMC_Geodatabase"), recursive = TRUE)
  file.remove(path_sgmc_zip)
  message("Done: SGMC geodatabase saved to data-raw")
} else {
  message("Already present, skipping: SGMC geodatabase")
}

# ---- LANID (irrigation) -----------------------------------------------------
# Archive bundles 2011-2017. Only lanid2017.tif is kept.
path_lanid_zip <- here("data-raw", "lanid2011-2017.zip")
path_lanid     <- here("data-raw", "lanid2017.tif")

if (!file.exists(path_lanid)) {
  message("Downloading LANID 2011-2017 archive (~345 MB) ...")
  download.file(
    url      = "https://zenodo.org/records/5003977/files/lanid2011-2017.zip",
    destfile = path_lanid_zip,
    mode     = "wb"
  )
  message("Extracting lanid2017.tif ...")
  unzip(path_lanid_zip, files = "lanid2017.tif", exdir = here("data-raw"))
  file.remove(path_lanid_zip)
  message("Done: LANID 2017 saved to data-raw")
} else {
  message("Already present, skipping: LANID 2017")
}

# ---- TWI (topographic wetness index) ----------------------------------------
# At ~41 GB the download eventually breaks,
# and curl -C - resumes a partial file. If complete, curl exits without action.
# If curl is unavailable, download manually from the DOI in the README.md and place 
# into data-raw.
path_twi <- here("data-raw", "CONUS_TWI_epsg5072_30m_unmasked.tif")

message("Checking / downloading CONUS TWI (~41 GB) ...")
exit_code <- system(paste(
  "curl -C - -L -o", shQuote(path_twi),
  shQuote(paste0(
    "https://zenodo.org/records/4460354/files/",
    "CONUS_TWI_epsg5072_30m_unmasked.tif?download=1"
  ))
))

# Exit 0 = ok. Exit 33 = server refused the resume range (file already complete).
if (!exit_code %in% c(0, 33)) {
  stop("curl exited with code ", exit_code,
       ". Re-run to resume, or download manually.")
}
message("Done: CONUS TWI present in data-raw/")