# ============================================================================
# compare_reference_years.R
# Compares two runs of the cwd_max pipeline that differ only in the length
# of the reference-years baseline used for ET gap-filling.
#
# Assumes combined multi-band exports (cwd_max + gap_filled_months) with
# the filenames below. Adjust paths/names as needed.
# ============================================================================

library(terra)
library(here)

# ── Load ─────────────────────────────────────────────────────────────────────
r_short <- rast(here("data", "cwd_test_subset_eel_ref3.tif"))   # 2020–2022
r_long  <- rast(here("data", "cwd_test_subset_eel_ref9.tif"))   # 2016–2024

cwd_s <- r_short[["cwd_max"]]
cwd_l <- r_long [["cwd_max"]]
gap_s <- r_short[["gap_filled_months"]]
gap_l <- r_long [["gap_filled_months"]]

# ── 1. NA coverage ───────────────────────────────────────────────────────────
cat("── NA coverage (cwd_max) ──\n")
n_total   <- ncell(cwd_s)
na_short  <- global(is.na(cwd_s), "sum")[1, 1]
na_long   <- global(is.na(cwd_l), "sum")[1, 1]
recovered <- global(is.na(cwd_s) & !is.na(cwd_l), "sum")[1, 1]
lost      <- global(!is.na(cwd_s) & is.na(cwd_l), "sum")[1, 1]

cat(sprintf("  Total cells:              %d\n", n_total))
cat(sprintf("  NA with 3-year reference: %d  (%.1f%%)\n", na_short, 100 * na_short / n_total))
cat(sprintf("  NA with 9-year reference: %d  (%.1f%%)\n", na_long,  100 * na_long  / n_total))
cat(sprintf("  Pixels recovered (NA→value): %d\n", recovered))
cat(sprintf("  Pixels lost     (value→NA):  %d\n\n", lost))

# ── 2. Value differences where both have values ──────────────────────────────
both_valid <- !is.na(cwd_s) & !is.na(cwd_l)
diff       <- cwd_l - cwd_s
diff_valid <- mask(diff, both_valid, maskvalues = c(0, NA))

cat("── Value differences (9yr − 3yr, mm, both valid) ──\n")
print(global(diff_valid, c("min", "max", "mean", "sd"), na.rm = TRUE))
cat("\n")

# ── 3. Differences conditional on gap-filling ────────────────────────────────
# At pixels that were gap-filled at least once (in either run), the
# reference-mean difference actually enters the CWD calculation — those
# are the pixels where the extra reference years should matter.
gap_any <- (gap_s > 0) | (gap_l > 0)
diff_gf <- mask(diff_valid, gap_any, maskvalues = c(0, NA))

cat("── Value differences AT gap-filled pixels ──\n")
print(global(diff_gf, c("min", "max", "mean", "sd"), na.rm = TRUE))
cat(sprintf("  Gap-filled pixels (either run): %d\n\n",
            global(gap_any, "sum", na.rm = TRUE)[1, 1]))

# ── 4. Plots ─────────────────────────────────────────────────────────────────
par(mfrow = c(2, 2))
plot(cwd_s, main = "cwd_max — 3yr reference",  colNA = "magenta")
plot(cwd_l, main = "cwd_max — 9yr reference",  colNA = "magenta")
plot(diff_valid, main = "Difference (9yr − 3yr)",
     col = hcl.colors(100, "Blue-Red"))
plot(is.na(cwd_s) & !is.na(cwd_l),
     main = "Pixels recovered by longer reference",
     col  = c("grey90", "darkgreen"))
par(mfrow = c(1, 1))