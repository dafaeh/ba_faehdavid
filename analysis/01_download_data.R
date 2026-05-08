# 01_download_data.R
# Downloads all external data required for the analysis.
# GEE exports (cwd rasters, elevation rasters) are excluded — these are
# downloaded manually from Google Drive and placed in data-raw/.
# Run once before any analysis script.
#
# Data source overview:
#
#   DOWNLOADED HERE:
#     - SGMC (geology): Horton et al. (2017), https://doi.org/10.5066/F7WH2N65
#     - LANID (irrigation): Xie et al. (2021), https://doi.org/10.5194/essd-13-5689-2021
#
#   STREAMED ON DEMAND (no setup required, accessed in analysis scripts):
#     - NLCD (land cover): Dewitz (2023), https://doi.org/10.5066/P9JZ7AO3
#       Accessed via /vsicurl/ from the FedData mirror on Google Cloud Storage.
#
#   DOWNLOADED MANUALLY (GEE exports, place in data-raw/):
#     - CWD rasters:       *_Xref.tif
#     - Elevation rasters: elevation_*.tif
#
#   HANDLED IN GEE (no local download required):
#     - FABDEM (terrain): Hawker et al. (2022), https://doi.org/10.1088/1748-9326/ac4d4f
#       TWI and elevation are computed in GEE and exported alongside CWD rasters.


library(here)

# Raise timeout for large file downloads (~400-500 MB).
options(timeout = 600)


# ---- Helper -----------------------------------------------------------------

# Skips download if the target file already exists. Prevents redundant
# re-downloads on repeated runs of this script.
download_if_missing <- function(path, url, description) {
  if (file.exists(path)) {
    message("Already present, skipping: ", description)
    return(invisible(NULL))
  }
  message("Downloading: ", description, " ...")
  download.file(url = url, destfile = path, mode = "wb")
  message("Done: ", description)
}


# ---- SGMC (geology) ---------------------------------------------------------
# Horton et al. (2017). The State Geologic Map Compilation (SGMC) geodatabase
# of the conterminous United States (ver. 1.1).
# https://doi.org/10.5066/F7WH2N65
#
# The file hash in the ScienceBase URL is content-based and stable across
# server changes. If the URL breaks, download manually from the DOI above
# and place the .gdb in data-raw/.
#
# Note: the ZIP is ~400 MB; download may take several minutes.

path_sgmc_zip <- here("data-raw", "USGS_SGMC_Geodatabase.zip")
path_sgmc_gdb <- here("data-raw", "USGS_StateGeologicMapCompilation_ver1.1.gdb")

if (!file.exists(path_sgmc_gdb)) {
  download_if_missing(
    path        = path_sgmc_zip,
    url         = paste0(
      "https://www.sciencebase.gov/catalog/file/get/",
      "5888bf4fe4b05ccb964bab9d",
      "?f=__disk__24/f6/e1/",
      "24f6e139c181fd9fe43df2aaf7f50b1c5b3b6297"
    ),
    description = "SGMC geodatabase (~400 MB)"
  )
  message("Unzipping SGMC geodatabase ...")
  unzip(path_sgmc_zip, exdir = here("data-raw"))
  file.rename(
    here("data-raw", "USGS_SGMC_Geodatabase",
         "USGS_StateGeologicMapCompilation_ver1.1.gdb"),
    path_sgmc_gdb
  )
  unlink(here("data-raw", "USGS_SGMC_Geodatabase"), recursive = TRUE)
  file.remove(path_sgmc_zip)
  message("Done.")
} else {
  message("Already present, skipping: SGMC geodatabase")
}


# ---- LANID (irrigation) -----------------------------------------------------
# Xie & Lark (2021). Landsat-based Irrigation Dataset (LANID).
# https://doi.org/10.5194/essd-13-5689-2021
# Data repository: https://doi.org/10.5281/zenodo.5003977
#
# The 2017 annual map is used as a proxy for irrigation status during the
# 2020-2022 CWD analysis period (most recent available year). Years 2011-2017
# are bundled in a single ZIP archive; only lanid2017.tif is extracted.
# The ZIP is deleted after extraction to save disk space (~361 MB).
#
# The global timeout set above is required for this ~345 MB file.

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
  message("Done: LANID 2017 saved to data-raw/")
} else {
  message("Already present, skipping: LANID 2017")
}


# ---- Done -------------------------------------------------------------------

message("\n01_download_data.R complete. All required files are in data-raw/.")
message("Reminder: GEE exports (CWD rasters, elevation rasters) must be")
message("placed in data-raw/ manually from Google Drive.")