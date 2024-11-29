#' @title plot.modelFits
#'
#' @description Plot function based on the ggplot2 package. Providing visualizations for each model and a average Fit.
#' Black lines show the fitted dose response models and an AIC based average model.
#' Dots indicate the posterior median and vertical lines show corresponding credible intervals (i.e. the variability of the posterior distribution of the respective dose group).
#' To assess the uncertainty of the model fit one can in addition visualize credible bands (default coloring as orange shaded areas).
#' The calculation of these bands is performed via the getBootstrapQuantiles() function.
#' The default setting is that these credible bands are not calculated.
#' @param x An object of type modelFits
#' @param gAIC Logical value indicating whether gAIC values are shown in the plot. Default TRUE
#' @param avg_fit Logical value indicating whether average fit is presented in the plot. Default TRUE
#' @param cr_intv Logical value indicating whether credible intervals are included in the plot. Default TRUE
#' @param alpha_CrI Numerical value of the width of the credible intervals. Default is set to 0.05 (i.e 95% CI are shown).
#' @param cr_bands Logical value indicating whether bootstrapped based credible bands are shown in the plot. Default FALSE
#' @param alpha_CrB Numerical vector of the width of the credible bands. Default is set to 0.05 and 0.5 (i.e 95% CB and median are shown).
#' @param n_bs_smpl Number of bootstrap samples being used. Default set to 1000.
#' @param acc_color Color of the credible bands. Default set to "orange"
#' @param ... optional parameter to be passed.
#' @examples
#' posterior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 1), sigma = 2),
#'                        DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 3, s = 1.2), sigma = 2),
#'                        DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 4, s = 1.5), sigma = 2) ,
#'                        DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 6, s = 1.2), sigma = 2) ,
#'                        DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 6.5, s = 1.1), sigma = 2))
#' models <- c("exponential", "linear")
#' dose_levels <- c(0, 1, 2, 4, 8)
#' fit <- getModelFits(models      = models,
#'                     posterior   = posterior_list,
#'                     dose_levels = dose_levels,
#'                     simple      = TRUE)
#'
#' plot(fit)
#' @return A ggplot2 object
#' @export
plot.modelFits <- function (

  x,
  gAIC      = TRUE,
  avg_fit   = TRUE,
  cr_intv   = TRUE,
  alpha_CrI = 0.05,
  cr_bands  = FALSE,
  alpha_CrB = c(0.05, 0.5),
  n_bs_smpl = 1e3,
  acc_color = "orange",
  ...

) {
  ## R CMD --as-cran appeasement
  .data <- NULL

  checkmate::check_class(x, "modelFits")
  checkmate::check_logical(gAIC)
  checkmate::check_logical(avg_fit)
  checkmate::check_logical(cr_intv)
  checkmate::check_double(alpha_CrI, lower = 0, upper = 1)
  checkmate::check_logical(cr_bands)
  checkmate::check_double(alpha_CrB, lower = 0, upper = 1, len = 2)
  checkmate::check_integer(n_bs_smpl, lower = 1, upper = Inf)
  checkmate::check_string(acc_color, na.ok = TRUE)

  plot_res   <- 1e2
  model_fits <- x

  dose_levels  <- model_fits[[1]]$dose_levels
  post_summary <- summary.postList(
    object = attr(model_fits, "posterior"),
    probs  = c(alpha_CrI / 2, 0.5, 1 - alpha_CrI / 2))
  dose_seq <- seq(from       = min(dose_levels),
                  to         = max(dose_levels),
                  length.out = plot_res)

  preds_models <- sapply(model_fits, predictModelFit, doses = dose_seq)
  model_names  <- names(model_fits)

  if (avg_fit) {

    mod_weights <- sapply(model_fits, function (x) x$model_weight)
    avg_preds   <- preds_models %*% mod_weights

    preds_models <- cbind(preds_models, avg_preds)
    model_names  <- c(model_names, "avgFit")

  }

  gg_data <- data.frame(
    doses  = rep(dose_seq, length(model_names)),
    fits   = as.vector(preds_models),
    models = rep(factor(model_names, levels = model_names),
                 each = plot_res))

  if (gAIC) {

    g_AICs     <- sapply(model_fits, function (x) x$gAIC)
    label_gAUC <- paste("AIC:", round(g_AICs, digits = 1))

    if (avg_fit) {

      mod_weights <- sort(mod_weights, decreasing = TRUE)
      paste_names <- names(mod_weights) |>
        gsub("exponential", "exp", x = _) |>
        gsub("quadratic",  "quad", x = _) |>
        gsub("linear",      "lin", x = _) |>
        gsub("logistic",    "log", x = _) |>
        gsub("sigEmax",    "sigE", x = _) |>
        gsub("betaMod",   "betaM", x = _) |>
        gsub("quadratic",  "quad", x = _)

      label_avg  <- paste0(paste_names, "=", round(mod_weights, 1),
                           collapse = ", ")
      label_gAUC <- c(label_gAUC, label_avg)

    }

  }

  if (cr_bands) {

    crB_data <- getBootstrapQuantiles(
      model_fits = model_fits,
      n_samples  = n_bs_smpl,
      quantiles  = c(0.5, sort(unique(c(alpha_CrB / 2, 1 - alpha_CrB / 2)))),
      avg_fit    = avg_fit,
      doses      = dose_seq)

    getInx <- function (alpha_CrB) {
      n        <- length(alpha_CrB)
      inx_list <- lapply(seq_len(n), function (i) c(i, 2 * n - i + 1) + 3)
      return (inx_list)}

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
        size = 3)

    }
    }

  if (cr_bands) {

    ## Bootstrapped Credible Bands
    for (inx in getInx(alpha_CrB)) {

      loop_txt <- paste0(
        "ggplot2::geom_ribbon(
          data    = crB_data,
          mapping = ggplot2::aes(x    = doses,
                                 ymin = crB_data[, ", inx[1], "],
                                 ymax = crB_data[, ", inx[2], "]),
          fill    = acc_color,
          alpha   = 0.25)")

      plts <- plts + eval(parse(text = loop_txt))

    }
    rm(inx)

    ## Bootstrapped Fit
    plts <- plts +
      ggplot2::geom_line(
        data    = crB_data,
        mapping = ggplot2::aes(.data$doses, .data$`50%`),
        color   = acc_color)

  }

  plts <- plts +
    ## Posterior Credible Intervals
    {if (cr_intv) {

      ggplot2::geom_errorbar(
        data    = data.frame(x    = dose_levels,
                             ymin = post_summary[, 3],
                             ymax = post_summary[, 5]),
        mapping = ggplot2::aes(x    = .data$x,
                               ymin = .data$ymin,
                               ymax = .data$ymax),
        width   = 0,
        alpha   = 0.5)

    }
    } +
    ## Posterior Medians
    ggplot2::geom_point(
      data    = data.frame(
        dose_levels = dose_levels,
        medians     = post_summary[, 4]),
      mapping = ggplot2::aes(.data$dose_levels, .data$medians),
      size    = 2) +
    ## Fitted Models
    ggplot2::geom_line(
      data    = gg_data,
      mapping = ggplot2::aes(.data$doses, .data$fits)) +
    ## Faceting
    ggplot2::facet_wrap(~ models)

  return (plts)

}
