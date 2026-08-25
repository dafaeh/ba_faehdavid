# Save all relevant plots
library(terra)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)


# ── Setup ───────────────────────────────────────────────────────────────────

# create figures directory if it doesn't exist
dir.create(here("figures"), showWarnings = FALSE)

# load data
eel_3ref     <- rast(here("data", "eel_subset_3ref_combined.tif"))
eel_9ref     <- rast(here("data", "eel_subset_9ref_combined.tif"))
app_3ref     <- rast(here("data", "appalachia_subset_3ref_snow.tif"))
app_9ref     <- rast(here("data", "appalachia_subset_9ref_snow.tif"))
app_no_snow  <- rast(here("data", "appalachia_subset_9ref_nosnow.tif"))


# ── Eel-Subset cwd, 3ref ────────────────────────────────────────────────────

png(here("figures", "eel_subset_cwd_3ref.png"), width = 800, height = 600)
plot(eel_3ref[["cwd_max"]],
     main  = "Eel-Subset: cwd_max (mm), 3 Referenzjahre",
     colNA = "magenta")
dev.off()


# ── Eel-Subset gap-filling, 3ref ────────────────────────────────────────────

png(here("figures", "eel_subset_gap_3ref.png"), width = 800, height = 600)
plot(eel_3ref[["gap_filled_months"]],
     main  = "Eel-Subset: gap_filled_months, 3 Referenzjahre",
     colNA = "magenta")
dev.off()


# ── Eel-Subset cwd, 9ref ────────────────────────────────────────────────────

png(here("figures", "eel_subset_cwd_9ref.png"), width = 800, height = 600)
plot(eel_9ref[["cwd_max"]],
     main  = "Eel-Subset: cwd_max (mm), 9 Referenzjahre",
     colNA = "magenta")
dev.off()


# ── Eel-Subset gap-filling, 9ref ────────────────────────────────────────────

png(here("figures", "eel_subset_gap_9ref.png"), width = 800, height = 600)
plot(eel_9ref[["gap_filled_months"]],
     main  = "Eel-Subset: gap_filled_months, 9 Referenzjahre",
     colNA = "magenta")
dev.off()


# ── Appalachen-Subset cwd, 3ref ─────────────────────────────────────────────

png(here("figures", "appalachia_subset_cwd_3ref.png"), width = 800, height = 600)
plot(app_3ref[["cwd_max"]],
     main  = "Appalachen-Subset: cwd_max (mm), 3 Referenzjahre",
     colNA = "magenta")
dev.off()


# ── Appalachen-Subset gap-filling, 3ref ─────────────────────────────────────

png(here("figures", "appalachia_subset_gap_3ref.png"), width = 800, height = 600)
plot(app_3ref[["gap_filled_months"]],
     main  = "Appalachen-Subset: gap_filled_months, 3 Referenzjahre",
     colNA = "magenta")
dev.off()


# ── Appalachen-Subset cwd, 9ref ─────────────────────────────────────────────

pdf(here("figures", "appalachia_subset_cwd_9ref.pdf"),
    width = 8, height = 6)
plot(app_9ref[["cwd_max"]],
     main  = "Appalachen-Subset: cwd_max (mm), 9 Referenzjahre",
     colNA = "magenta")
dev.off()


# ── Appalachen-Subset gap-filling, 9ref ─────────────────────────────────────

png(here("figures", "appalachia_subset_gap_9ref.png"), width = 800, height = 600)
plot(app_9ref[["gap_filled_months"]],
     main  = "Appalachen-Subset: gap_filled_months, 9 Referenzjahre",
     colNA = "magenta")
dev.off()


# ── Schneeeffekt-Differenzkarte ─────────────────────────────────────────────

diff <- app_9ref[["cwd_max"]] - app_no_snow[["cwd_max"]]

abs_max <- max(abs(global(diff, c("min", "max"), na.rm = TRUE)[1, ]))

png(here("figures", "appalachia_subset_snow_effect_map.png"),
    width = 800, height = 600)
plot(diff,
     main  = "Appalachen-Subset: Effekt des Schneemodells auf cwd_max (mm)",
     col   = hcl.colors(100, "Blue-Red"),
     range = c(-abs_max, abs_max))
dev.off()


# ── Schneeeffekt nach Höhenstufe ────────────────────────────────────────────

srtm <- rast("/vsicurl/https://opentopography.s3.sdsc.edu/raster/SRTM_GL1/SRTM_GL1_srtm.vrt")
elev <- crop(srtm, project(ext(diff), from = crs(diff), to = crs(srtm)))
elev <- project(elev, diff, method = "bilinear")

samp <- spatSample(c(diff, elev), size = 50000, method = "regular",
                   na.rm = TRUE, as.df = TRUE)
names(samp) <- c("diff", "elev")

samp <- samp |>
  drop_na() |>
  mutate(bin = cut(elev,
                   breaks = seq(0, 1500, 100),
                   labels = paste0(seq(0, 1400, 100), "-", seq(100, 1500, 100))))

agg <- samp |>
  group_by(bin) |>
  summarise(mean_diff = mean(diff))

p <- ggplot(agg, aes(x = bin, y = mean_diff)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0) +
  labs(
    title = "Appalachen-Subset — Schneeeffekt nach Höhenstufe",
    x     = "Höhe (m)",
    y     = "mittlere Differenz (mm)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here("figures", "appalachia_subset_snow_effect_by_elevation.png"),
       plot = p, width = 8, height = 5, dpi = 150)