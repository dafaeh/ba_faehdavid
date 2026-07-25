# save_fig.R
# Writes a ggplot to fig/ as both PDF and PNG under the same base name.

save_fig <- function(plot, name, width = 16, height = 12, dpi = 600) {
  ggplot2::ggsave(
    filename = here::here("fig", paste0(name, ".pdf")),
    plot     = plot,
    width    = width,
    height   = height,
    units    = "cm"
  )
  ggplot2::ggsave(
    filename = here::here("fig", paste0(name, ".png")),
    plot     = plot,
    width    = width,
    height   = height,
    units    = "cm",
    dpi      = dpi
  )
  invisible(plot)
}