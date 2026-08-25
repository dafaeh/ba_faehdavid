library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(here)

set.seed(42)

FORMATION_MAP <- list(
  "Edwards Ls."   = c("Edwards Limestone",
                      "Edwards and Comanche Peak Limestones, undivided"),
  "Glen Rose Ls." = "Glen Rose Limestone"
)

pal <- c("Edwards Ls." = "#d73027", "Glen Rose Ls." = "#fc8d59")
FORMATION_LEVELS <- names(FORMATION_MAP)
all_target_names <- unlist(FORMATION_MAP, use.names = FALSE)

# ---- Raster -----------------------------------------------------------------

r     <- terra::rast(here("data", "edwards_plateau_9ref.tif"))
r_cwd <- r[["cwd_max"]]

# ---- SGMC -------------------------------------------------------------------

gdb_path  <- here("data-raw", "USGS_StateGeologicMapCompilation_ver1.1.gdb")
aoi_wgs84 <- sf::st_as_sfc(sf::st_bbox(
  c(xmin = -100.27776, ymin = 28.78932,
    xmax =  -96.21282, ymax = 32.66018),
  crs = sf::st_crs(4326)
))
gdb_crs <- sf::st_crs(
  sf::st_read(gdb_path, query = "SELECT * FROM SGMC_Geology LIMIT 1", quiet = TRUE)
)

sgmc <- sf::st_read(
  dsn        = gdb_path,
  layer      = "SGMC_Geology",
  wkt_filter = sf::st_as_text(sf::st_transform(aoi_wgs84, gdb_crs)),
  quiet      = TRUE
) |>
  sf::st_transform(crs = 5070) |>
  dplyr::filter(UNIT_NAME %in% all_target_names) |>
  dplyr::mutate(
    formation = dplyr::case_when(
      UNIT_NAME %in% FORMATION_MAP[["Edwards Ls."]]   ~ "Edwards Ls.",
      UNIT_NAME %in% FORMATION_MAP[["Glen Rose Ls."]] ~ "Glen Rose Ls.",
      UNIT_NAME %in% FORMATION_MAP[["Hensell Sand"]]  ~ "Hensell Sand"
    ),
    formation = factor(formation, levels = FORMATION_LEVELS)
  )

# ---- Rasterize --------------------------------------------------------------

sgmc_vect          <- terra::vect(sgmc)
sgmc_vect$form_int <- as.integer(sgmc$formation)
form_rast          <- terra::rasterize(sgmc_vect, r_cwd, field = "form_int")

# ---- Sample formation raster, then extract CWD at points -------------------

n_sample <- 50000

sample_pts <- terra::spatSample(
  form_rast,
  size   = n_sample,
  method = "stratified",
  na.rm  = TRUE,
  xy     = TRUE,
  as.df  = TRUE
)

cwd_at_pts <- terra::extract(r, sample_pts[, c("x", "y")])

df <- dplyr::tibble(
  formation         = sample_pts$form_int,
  cwd_max           = cwd_at_pts$cwd_max,
  gap_filled_months = cwd_at_pts$gap_filled_months
) |>
  dplyr::mutate(
    formation = factor(
      formation,
      levels = seq_along(FORMATION_LEVELS),
      labels = FORMATION_LEVELS
    )
  ) |>
  dplyr::filter(!is.na(cwd_max))

message("Valid pixels by formation:")
print(dplyr::count(df, formation))

# ---- Konsole ----------------------------------------------------------------

df |>
  dplyr::group_by(formation) |>
  dplyr::summarise(
    n      = dplyr::n(),
    median = median(cwd_max),
    mean   = mean(cwd_max),
    sd     = sd(cwd_max),
    .groups = "drop"
  ) |>
  print(width = Inf)

# ---- Plots ------------------------------------------------------------------

pal <- c("Edwards Ls." = "#d73027", "Glen Rose Ls." = "#fc8d59",
         "Hensell Sand" = "#1a9850")

print(
  ggplot(df, aes(x = formation, y = cwd_max, fill = formation)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
    scale_fill_manual(values = pal) +
    labs(x = NULL, y = expression(paste(CWD[max], " (mm)"))) +
    theme_classic() +
    theme(legend.position = "none")
)

print(
  ggplot(df, aes(x = cwd_max, fill = formation, colour = formation)) +
    geom_density(alpha = 0.4) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_colour_manual(values = pal, name = NULL) +
    labs(x = expression(paste(CWD[max], " (mm)")), y = "Density") +
    theme_classic() +
    theme(legend.position = c(0.78, 0.82))
)