# ============================================================================
# compare_snow_vs_no_snow.R
# Compares two runs of the cwd_max pipeline over the Appalachian subset
# that differ only in TEMP_THRESHOLD:
#   - snow run:    TEMP_THRESHOLD =  1.0  (snow model active)
#   - no-snow run: TEMP_THRESHOLD = -999  (snow model effectively disabled,
#                                          all precipitation = rain)
#
# Goal: isolate the effect of the snow model on cwd_max.
# Hypothesis: where snow accumulates over winter and melts in spring, the
# snow run should show HIGHER cwd_max than the no-snow run, because winter
# precipitation is delayed and not available to offset summer ET.
# ============================================================================

library(terra)
library(here)

# ── Load ─────────────────────────────────────────────────────────────────────
r_snow    <- rast(here("data", "explorative", "appalachia_subset_9ref_snow.tif"))
r_no_snow <- rast(here("data", "explorative", "appalachia_subset_9ref_nosnow.tif"))

cwd_snow    <- r_snow   [["cwd_max"]]
cwd_no_snow <- r_no_snow[["cwd_max"]]

# ── Geometry sanity ──────────────────────────────────────────────────────────
stopifnot(ext(cwd_snow) == ext(cwd_no_snow),
          all(res(cwd_snow) == res(cwd_no_snow)),
          same.crs(cwd_snow, cwd_no_snow))

# ── Difference: snow − no_snow ───────────────────────────────────────────────
# Positive = snow model INCREASES cwd_max (winter precip held back)
# Negative = snow model DECREASES cwd_max (would be unexpected)
diff <- cwd_snow - cwd_no_snow

cat("── Summary of cwd_max (snow run) ──\n")
print(global(cwd_snow, c("min", "max", "mean", "sd"), na.rm = TRUE))
cat("\n── Summary of cwd_max (no-snow run) ──\n")
print(global(cwd_no_snow, c("min", "max", "mean", "sd"), na.rm = TRUE))

cat("\n── Difference (snow − no_snow) ──\n")
print(global(diff, c("min", "max", "mean", "sd"), na.rm = TRUE))

# Distribution of the snow-model effect
n_total       <- global(!is.na(diff), "sum")[1, 1]
n_unchanged   <- global(abs(diff) < 1e-3, "sum", na.rm = TRUE)[1, 1]
n_higher_snow <- global(diff >  1e-3,     "sum", na.rm = TRUE)[1, 1]
n_lower_snow  <- global(diff < -1e-3,     "sum", na.rm = TRUE)[1, 1]

cat(sprintf("\nValid pixels:                   %d\n", n_total))
cat(sprintf("Effectively unchanged (|d|<1mm): %d  (%.1f%%)\n",
            n_unchanged,   100 * n_unchanged   / n_total))
cat(sprintf("Snow run HIGHER  (d > 1 mm):     %d  (%.1f%%)\n",
            n_higher_snow, 100 * n_higher_snow / n_total))
cat(sprintf("Snow run LOWER   (d < -1 mm):    %d  (%.1f%%)\n",
            n_lower_snow,  100 * n_lower_snow  / n_total))

# ── Plots ────────────────────────────────────────────────────────────────────
par(mfrow = c(2, 2))

plot(cwd_snow,    main = "cwd_max — with snow model",  colNA = "magenta")
plot(cwd_no_snow, main = "cwd_max — no snow",          colNA = "magenta")

# Symmetric diverging scale around zero so the colour map is interpretable.
abs_max <- max(abs(global(diff, c("min", "max"), na.rm = TRUE)[1, ]))
plot(diff,
     main = "Difference (snow − no_snow), mm",
     col  = hcl.colors(100, "Blue-Red"),
     range = c(-abs_max, abs_max))

# Histogram of the differences for a quick look at the distribution.
hist(values(diff),
     breaks = 80,
     main   = "Distribution of differences",
     xlab   = "cwd_max(snow) − cwd_max(no_snow)  [mm]")
abline(v = 0, lty = 2)

par(mfrow = c(1, 1))