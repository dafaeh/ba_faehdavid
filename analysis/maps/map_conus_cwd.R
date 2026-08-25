# ---- Setup ------------------------------------------------------------------
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyterra)
library(scico)
library(scales)
library(here)

# ---- Load data --------------------------------------------------------------
r_conus      <- terra::rast(here("data", "cwd_conus.tif"))
r_cwd        <- r_conus[["cwd_max"]]
r_gap_filled <- r_conus[["gap_filled_months"]]

# ---- Downsample for plotting ------------------------------------------------
# 30 m -> 300 m (factor 10). Continuous variables -> mean aggregation.
r_cwd_agg        <- terra::aggregate(r_cwd,        fact = 10, fun = "mean")
r_gap_filled_agg <- terra::aggregate(r_gap_filled, fact = 10, fun = "mean")

# ---- Map: CWDmax ------------------------------------------------------------
p_map_cwd <- ggplot() +
  tidyterra::geom_spatraster(data = r_cwd_agg) +
  scale_fill_viridis_c(
    name      = expression(CWD[max]~"[mm]"),
    option    = "viridis",
    direction = -1,              # low = yellow, high = blue/purple
    limits    = c(0, 840),
    oob       = scales::squish,  # values > 840 shown as darkest colour
    na.value  = NA
  ) +
  labs(
    title = expression("CWD"[max]~"across the contiguous United States [mm]")
  ) +
  theme_classic()
p_map_cwd

# ---- Map: gap-filled months -------------------------------------------------
# Separate palette family from the CWDmax map so the two are not confused. 
# Range 0-36: sum of gap-filled months per pixel across the
# three analysis years 
p_map_gap <- ggplot() +
  tidyterra::geom_spatraster(data = r_gap_filled_agg) +
  scico::scale_fill_scico(
    name      = "Gap-filled\nmonths [count]",
    palette   = "lajolla",
    direction = -1,              # low = light, high = dark red/brown
    limits    = c(0, 36),
    oob       = scales::squish,
    na.value  = NA
  ) +
  labs(
    title = "Gap-filled months across the contiguous United States"
  ) +
  theme_classic()
p_map_gap

# ---- Sample for the distribution plot ---------------------------------------
# na.rm = FALSE plus drop_na() afterwards: na.rm = TRUE makes spatSample
# oversample repeatedly until enough valid pixels are found, which is
# prohibitively slow on a raster of this size.
set.seed(42)
df_sample <- terra::spatSample(
  r_conus,
  size   = 100000,
  method = "random",
  na.rm  = FALSE,
  as.df  = TRUE
) |>
  tibble::as_tibble() |>
  tidyr::drop_na()

# ---- Distribution: gap-filled months ----------------------------------------
# Discrete count (0-36) -> one bar per integer value, relative frequency.
p_bar_gap <- ggplot(df_sample, aes(x = gap_filled_months)) +
  geom_bar(
    aes(y = after_stat(count / sum(count))),
    fill = "grey30", width = 0.8
  ) +
  scale_x_continuous(breaks = seq(0, 36, 6)) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Gap-filled months (count)",
    y = "Share of pixels"
  ) +
  theme_classic()
p_bar_gap

# ---- Combined figure: gap-filled months -------------------------------------
# Map and distribution as one thesis figure. Titles are dropped here, the
# caption describes both panels.
p_gap_combined <- cowplot::plot_grid(
  p_map_gap + labs(title = NULL),
  p_bar_gap,
  ncol        = 1,
  labels      = "auto",
  align       = "v",
  axis        = "lr",          # panel edges stay flush despite the map legend
  rel_heights = c(1.6, 1)
)
p_gap_combined

# ---- Save -------------------------------------------------------------------
ggsave(here("fig", "map_cwd_max_conus.pdf"), p_map_cwd,
       width = 20, height = 14, units = "cm")
ggsave(here("fig", "map_cwd_max_conus.png"), p_map_cwd,
       width = 20, height = 14, units = "cm", dpi = 600)

ggsave(here("fig", "map_gap_filled_conus.pdf"), p_map_gap,
       width = 20, height = 14, units = "cm")
ggsave(here("fig", "map_gap_filled_conus.png"), p_map_gap,
       width = 20, height = 14, units = "cm", dpi = 600)

ggsave(here("fig", "bar_gap_filled_conus.pdf"), p_bar_gap,
       width = 16, height = 10, units = "cm")
ggsave(here("fig", "bar_gap_filled_conus.png"), p_bar_gap,
       width = 16, height = 10, units = "cm", dpi = 600)

ggsave(here("fig", "gap_filled_conus_combined.pdf"), p_gap_combined,
       width = 18, height = 16, units = "cm")
ggsave(here("fig", "gap_filled_conus_combined.png"), p_gap_combined,
       width = 18, height = 16, units = "cm", dpi = 600)