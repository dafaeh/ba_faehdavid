# =============================================================================
# correlate_twi_cwd_eel_subset.R
#
# Quick exploratory correlation analysis between Hoylman TWI and the
# project's cwd_max raster, both for the Eel River test subset.
#
# Paths are resolved relative to the project root via here::here(), so
# this script runs regardless of the current working directory as long
# as it lives somewhere inside the project.
# =============================================================================

library(terra)
library(here)

# ── Config ───────────────────────────────────────────────────────────────────
TWI_PATH <- here("data", "twi_eel_subset_hoylman_5070_30m.tif")
CWD_PATH <- here("figures", "eel_subset_cwd_3ref.tif")

# Maximum allowed number of gap-filled months across the 3-year analysis
# period (max possible = 36). Set to NA to disable the filter.
MAX_GAP_FILLED_MONTHS <- 6

# ── Load ─────────────────────────────────────────────────────────────────────
stopifnot(
  "TWI raster not found"     = file.exists(TWI_PATH),
  "cwd_max raster not found" = file.exists(CWD_PATH)
)

twi <- rast(TWI_PATH)
cwd <- rast(CWD_PATH)

cat("\nTWI raster:\n");    print(twi)
cat("\ncwd_max raster:\n"); print(cwd)

# Extract individual bands.
cwd_max <- cwd[["cwd_max"]]
gap     <- cwd[["gap_filled_months"]]

# ── Align TWI onto the cwd_max grid ──────────────────────────────────────────
# Both should already be EPSG:5070 30 m, but the grid origins are likely
# offset because they came from independent reprojection chains. resample()
# brings TWI onto the exact cwd_max grid using bilinear interpolation
# (TWI is continuous).
if (!compareGeom(twi, cwd_max, stopOnError = FALSE)) {
  cat("\nGrids differ — resampling TWI onto the cwd_max grid.\n")
  twi <- resample(twi, cwd_max, method = "bilinear")
}

# ── Quality mask ─────────────────────────────────────────────────────────────
if (!is.na(MAX_GAP_FILLED_MONTHS)) {
  qmask   <- gap <= MAX_GAP_FILLED_MONTHS
  cwd_max <- mask(cwd_max, qmask, maskvalue = FALSE)
  cat(sprintf(
    "\nQuality mask applied: keeping pixels with <= %d gap-filled months.\n",
    MAX_GAP_FILLED_MONTHS
  ))
}

# ── Pull values into a data.frame ────────────────────────────────────────────
df <- as.data.frame(c(twi, cwd_max), na.rm = TRUE)
names(df) <- c("twi", "cwd_max")
cat(sprintf("\nValid pixel pairs: %d\n", nrow(df)))

# ── Correlation ──────────────────────────────────────────────────────────────
r_pearson  <- cor(df$twi, df$cwd_max, method = "pearson")
r_spearman <- cor(df$twi, df$cwd_max, method = "spearman")

cat(sprintf("\nPearson  r = %+.3f\n", r_pearson))
cat(sprintf("Spearman r = %+.3f\n",   r_spearman))

# Hypothesis H2 (valley bottoms accumulate higher CWD via lateral subsurface
# flow): would predict a *positive* correlation. A negative or near-zero
# correlation in the test subset is informative — possibly because shading
# effects or vegetation differences dominate in this small AOI.

# ── Diagnostic plots ─────────────────────────────────────────────────────────
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Density scatter — clearer than overplotting at this pixel count.
smoothScatter(
  df$twi, df$cwd_max,
  xlab    = "TWI",
  ylab    = "CWDmax (mm)",
  main    = sprintf("Pearson %.2f  /  Spearman %.2f",
                    r_pearson, r_spearman),
  colramp = colorRampPalette(c("white", "steelblue", "black"))
)
abline(lm(cwd_max ~ twi, data = df), col = "firebrick", lwd = 1.5)

# Side-by-side maps (rescaled — cwd_max only after the quality mask).
plot(c(twi, cwd_max),
     main = c("TWI", "cwd_max (masked)"),
     col  = hcl.colors(100, "Viridis"))

par(op)