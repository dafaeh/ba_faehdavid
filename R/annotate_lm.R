# R/annotate_lm.R
# Adds a label with R-squared, n, and (optionally) the linear slope to a
# ggplot scatter/hexbin plot that already contains an lm trend line.

annotate_lm <- function(p, model, slope_unit = NULL, digits_r2 = 2,
                        digits_slope = 2, size = 3.2) {
  glance_mod <- broom::glance(model)
  tidy_mod   <- broom::tidy(model)
  
  r2 <- round(glance_mod$r.squared, digits_r2)
  n  <- format(glance_mod$nobs, big.mark = "'")
  
  label <- paste0("R\u00b2 = ", r2, ", n = ", n)
  
  # Slope is optional
  if (!is.null(slope_unit)) {
    slope <- round(tidy_mod$estimate[2], digits_slope)
    label <- paste0(label, "\nslope = ", slope, " ", slope_unit)
  }
  
  p + annotate(
    "label",
    x = Inf, y = Inf, hjust = 1.05, vjust = 1.1,
    label = label, size = size,
    label.size    = 0.3,    
    label.padding = unit(0.4, "lines"), 
    fill = scales::alpha("white", 0.9)  
  )
}