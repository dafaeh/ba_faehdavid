# Quick raster inspection
library(terra)
library(here)


# ── Multi-band file ─────────────────────────────────────────────────────────

# load data
file <- "cwd_r01_c02-0000000000-0000023296.tif"
band <- "cwd_max"          # or "gap_filled_months"

r <- rast(here("data", file))[[band]]

# basic info
print(r)

# NA statistics
n_total <- ncell(r)
n_na    <- global(is.na(r), "sum")[1, 1]
n_valid <- n_total - n_na

cat("Total cells:   ", n_total, "\n")
cat("Valid cells:   ", n_valid, "\n")
cat("NA cells:      ", n_na, "\n")
cat("NA percentage: ", round(100 * n_na / n_total, 3), "%\n")

# value range
print(global(r, c("min", "max", "mean"), na.rm = TRUE))

# plot raster
plot(r, main = paste(band, "—", file), colNA = "white", range = c(0, 1400))

# plot NA mask
#plot(is.na(r), main = "NA mask (1 = NA, 0 = valid)")


# ── Single-band file ─────────────────────────────────────────────────

# load data
file <- "twi_appalachia.tif"

r <- rast(here("data", file))

# basic info
print(r)

# NA statistics
n_total <- ncell(r)
n_na    <- global(is.na(r), "sum")[1, 1]
n_valid <- n_total - n_na

cat("Total cells:   ", n_total, "\n")
cat("Valid cells:   ", n_valid, "\n")
cat("NA cells:      ", n_na, "\n")
cat("NA percentage: ", round(100 * n_na / n_total, 3), "%\n")

# value range
print(global(r, c("min", "max", "mean"), na.rm = TRUE))

# plot raster
plot(r, main = file, colNA = "white")

# plot NA mask
plot(is.na(r), main = "NA mask (1 = NA, 0 = valid)")
