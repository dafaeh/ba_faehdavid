# lanid_filter.R
# Script for excluding irrigated pixels from an analysis:
# a pixel is kept only if LANID marks it as non-irrigated.

# mask_irrigated() for scripts that filter on the
# raster before sampling, drop_irrigated() for scripts that sample first and
# filter with LANID afterwards.

# Raster version
mask_irrigated <- function(x, lanid) {
  keep <- terra::ifel(is.na(lanid), 0, lanid != 1)
  terra::mask(x, keep, maskvalue = 0)
}

# Data frame version
drop_irrigated <- function(df, col = "lanid") {
  dplyr::filter(df, !is.na(.data[[col]]) & .data[[col]] != 1)
}