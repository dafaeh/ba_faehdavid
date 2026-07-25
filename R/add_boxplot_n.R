# R/add_boxplot_n.R
# Computes n per group and returns a geom_text layer for annotating
# boxplots.


add_boxplot_n <- function(df, x_var, y_var, y_pos = NULL, size = 3) {
  df_n <- df |>
    dplyr::group_by(.data[[x_var]]) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  if (is.null(y_pos)) {
    y_pos <- max(df[[y_var]], na.rm = TRUE) * 1.05
  }
  
  geom_text(
    data        = df_n,
    mapping     = ggplot2::aes(x = .data[[x_var]], y = y_pos,
                               label = paste0("n = ", n)),
    inherit.aes = FALSE,
    size        = size
  )
}