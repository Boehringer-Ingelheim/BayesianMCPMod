plot.modelFits <- function (
    
  model_fits,
  CrI     = FALSE,
  gAIC    = TRUE,
  avg_fit = TRUE
  
) {
  
  plot_resolution <- 1e3
  
  dose_levels  <- model_fits[[1]]$dose_levels
  post_summary <- summary.postList(attr(model_fits, "posterior"))
  doses        <- seq(from = min(dose_levels),
                      to   = max(dose_levels), length.out = plot_resolution)
  
  preds_models <- sapply(model_fits, predictModelFit, doses = doses)
  model_names  <- names(model_fits)
  
  if (avg_fit) {
    
    mod_weigts <- sapply(model_fits, function (x) x$model_weight)
    avg_mod    <- preds_models %*% mod_weigts
    
    preds_models <- cbind(preds_models, avg_mod)
    model_names  <- c(model_names, "averageModel")
    
  }
  
  gg_data      <- data.frame(
    dose_levels = rep(doses, length(model_names)),
    fits        = as.vector(preds_models),
    models      = rep(factor(model_names,
                             levels = c("linear", "emax", "exponential",
                                        "sigEmax", "logistic", "quadratic",
                                        "averageModel")),
                      each = plot_resolution))
  
  if (gAIC) {
    
    g_AICs     <- sapply(model_fits, function (x) x$gAIC)
    label_gAUC <- paste("AIC:", round(g_AICs, digits = 1))
    
    if (avg_fit) {
      
      mod_weigts  <- sort(mod_weigts, decreasing = TRUE)
      paste_names <- names(mod_weigts) |>
        gsub("exponential", "exp", x = _) |>
        gsub("quadratic",  "quad", x = _) |>
        gsub("linear",      "lin", x = _) |>
        gsub("logistic",    "log", x = _) |>
        gsub("sigEmax",    "sigE", x = _)
      
      label_avg  <- paste0(paste_names, "=", round(mod_weigts, 1),
                           collapse = ", ")
      label_gAUC <- c(label_gAUC, label_avg)
      
    }
    
  }
  
  plts <- ggplot2::ggplot() + 
    ## Layout etc.
    ggplot2::theme_bw() + 
    ggplot2::labs(x = "Dose",
                  y = "Model Fits") +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ## gAIC
    {if (gAIC) {
      ggplot2::geom_text(
        data    = data.frame(
          models = unique(gg_data$models),
          label  = label_gAUC),
        mapping = ggplot2::aes(label = label_gAUC),
        x = -Inf, y = Inf, hjust = "inward", vjust = "inward",
        size = 3)}
    } + 
    ## Posterior Credible Intervals
    {if (CrI) {
      ggplot2::geom_errorbar(
        data    = data.frame(x    = dose_levels,
                             ymin = post_summary[, 3],
                             ymax = post_summary[, 5]),
        mapping = ggplot2::aes(x    = x,
                               ymin = ymin,
                               ymax = ymax),
        width = 0, alpha = 0.5)}
    } + 
    ## Posterior Medians
    ggplot2::geom_point(
      data    = data.frame(dose_levels = dose_levels,
                           fits        = post_summary[, 4]),
      mapping = ggplot2::aes(dose_levels, fits),
      size    = 2) +
    ## Fitted Models
    ggplot2::geom_line(
      data    = gg_data,
      mapping = ggplot2::aes(dose_levels, fits)) + 
    ## Faceting
    ggplot2::facet_wrap(~ models)
  
  return (plts)
  
}