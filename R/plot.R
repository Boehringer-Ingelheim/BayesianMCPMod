#' @title plot.modelFits
#'
#' @description Plot function based on the ggplot2 package. Providing visualizations for each model and a average Fit.
#' Black lines show the fitted dose response models and an AIC based average model.
#' Dots indicate the posterior median and vertical lines show corresponding credible intervals (i.e. the variability of the posterior distribution of the respective dose group).
#' To assess the uncertainty of the model fit one can in addition visualize credible bands (default coloring as orange shaded areas).
#' The calculation of these bands is performed via the getBootstrapQuantiles() function.
#' The default setting is that these credible bands are not calculated.
#' @param x An object of type modelFits
#' @param probability_scale A boolean to specify if the trial has a continuous or a binary outcome. Setting to TRUE will transform the output from the logit scale to the probability scale, which can be desirable for a binary outcome. Default `attr(x, "probability_scale")`.
#' @param gAIC Logical value indicating whether gAIC values are shown in the plot. Default TRUE
#' @param cr_intv Logical value indicating whether credible intervals are included in the plot. Default TRUE
#' @param alpha_CrI Numerical value of the width of the credible intervals. Default is set to 0.05 (i.e 95% CI are shown).
#' @param cr_bands Logical value indicating whether bootstrapped based credible bands are shown in the plot. Default FALSE
#' @param alpha_CrB Numerical vector of the width of the credible bands. Default is set to 0.05 and 0.5 (i.e 95% CB and 50% CB  are shown).
#' @param n_bs_smpl Number of bootstrap samples being used. Default 1000.
#' @param acc_color Color of the credible bands. Default "orange".
#' @param plot_res Number of plotted doses within the range of the dose levels, i.e., the resolution of the plot. Default 100.
#' @param plot_diffs Logical flag. Plot differences to the control group. Default FALSE.
#' @param ... optional parameter to be passed to plot().
#' @examples
#' posterior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 1), sigma = 2),
#'                        DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 3, s = 1.2), sigma = 2),
#'                        DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 4, s = 1.5), sigma = 2) ,
#'                        DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 6, s = 1.2), sigma = 2) ,
#'                        DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 6.5, s = 1.1), sigma = 2))
#' models <- c("exponential", "linear", "emax")
#' dose_levels <- c(0, 1, 2, 4, 8)
#' model_fits <- getModelFits(models      = models,
#'                            posterior   = posterior_list,
#'                            dose_levels = dose_levels,
#'                            simple       = TRUE)
#'
#' plot(model_fits)
#' 
#' # plot with credible bands
#' plot(model_fits,
#'      cr_bands          = TRUE,
#'      plot_diffs        = FALSE,
#'      probability_scale = FALSE,
#'      n_bs_smpl         = 1e2) # speeding up example run-time
#' 
#' @return A ggplot2 object
#' @export
plot.modelFits <- function (

  x,
  probability_scale = attr(x, "probability_scale"),
  gAIC              = TRUE,
  cr_intv           = TRUE,
  alpha_CrI         = 0.05,
  cr_bands          = FALSE,
  alpha_CrB         = c(0.05, 0.2),
  n_bs_smpl         = 1e3,
  acc_color         = "orange",
  plot_res          = 1e2,
  plot_diffs        = FALSE,
  ...

) {
  
  ## R CMD --as-cran appeasement
  .data <- q_prob <- q_val <- model <- sample_type <- NULL

  checkmate::assert_logical(gAIC)
  checkmate::assert_logical(cr_intv)
  checkmate::assert_double(alpha_CrI, lower = 0, upper = 1)
  checkmate::assert_logical(cr_bands)
  checkmate::assert_double(alpha_CrB, lower = 0, upper = 1, len = 2)
  checkmate::assert_integerish(n_bs_smpl, lower = 1, upper = Inf)
  checkmate::assert_string(acc_color, na.ok = TRUE)
  checkmate::assert_integerish(plot_res, lower = 1, upper = Inf)
  checkmate::assert_flag(probability_scale, null.ok = TRUE)
  checkmate::assert_flag(plot_diffs)
  
  if (is.null(probability_scale)) probability_scale <- FALSE

  model_fits  <- x
  model_names <- names(model_fits)
  avg_fit     <- "avgFit" %in% model_names

  dose_levels <- model_fits[[1]]$dose_levels
  dose_seq    <- seq(from       = min(dose_levels),
                     to         = max(dose_levels),
                     length.out = plot_res)
  
  preds_models <- do.call(c, lapply(
    predict.modelFits(model_fits, doses = dose_seq,
                      probability_scale = probability_scale),
    function (preds) if (plot_diffs) preds - preds[1] else preds))

  quantile_probs <- c(alpha_CrI / 2, 0.5, 1 - alpha_CrI / 2)
  if (plot_diffs) {
    
    if (probability_scale) {
      
      ## Sampling from posterior distributions
      post_samples   <- RBesT::inv_logit(sapply(attr(model_fits, "posterior"), RBesT::rmix, n = n_bs_smpl))
      post_quantiles <- apply(post_samples - post_samples[, 1], 2, stats::quantile, probs = quantile_probs)
      
    } else {
      
      post_quantiles <- sapply(attr(model_fits, "posterior"),
                               RBesT::qmixdiff,
                               mix2 = attr(model_fits, "posterior")[[1]],
                               p    = quantile_probs)
      post_quantiles[c(1, 3), 1] <- 0
      
    }
    
  } else {
    
    post_quantiles <- sapply(attr(model_fits, "posterior"),
                             RBesT::qmix,
                             p = quantile_probs)
    
    if (probability_scale) post_quantiles <- RBesT::inv_logit(post_quantiles)
    
  }
  
  ## ggplot data frama
  gg_data <- data.frame(
    doses  = rep(dose_seq, length(model_names)),
    fits   = as.vector(preds_models),
    models = rep(factor(model_names, levels = model_names),
                 each = plot_res))

  if (gAIC) {

    g_AICs     <- sapply(model_fits, function (x) x$gAIC)
    label_gAUC <- paste("AIC:", round(g_AICs, digits = 1))

    if (avg_fit) {
      
      mod_weights <- sapply(model_fits, function (y) y$model_weight)[-1]

      mod_weights <- sort(mod_weights, decreasing = TRUE)
      paste_names <- shortenModelNames(names(mod_weights))

      label_avg  <- paste0(paste_names, "=", round(mod_weights, 1),
                           collapse = ", ")
      label_gAUC <- paste0(" ", c(label_avg, label_gAUC[-1]))

    }

  }

  if (cr_bands) {
  
    quantile_pairs <- matrix(
      data = c(alpha_CrB / 2, rep(1, times = length(alpha_CrB))),
      ncol = 2)
    quantile_pairs[, 2] <- quantile_pairs[, 2] - quantile_pairs[, 1]

    crB_data <- getBootstrapQuantiles(
      model_fits        = model_fits,
      n_samples         = n_bs_smpl,
      quantiles         = sort(unique(c(0.5, as.vector(quantile_pairs)))),
      doses             = dose_seq,
      probability_scale = probability_scale) |>
      dplyr::filter(sample_type == dplyr::if_else(plot_diffs, "diff", "abs")) |>
      tidyr::pivot_wider(names_from = q_prob, values_from = q_val) |>
      dplyr::rename(models = model)

  }

  plts <- ggplot2::ggplot() +
    ## Layout etc.
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Dose",
                  y = "Model Fits") +
    # ggplot2::scale_x_continuous(breaks = c(dose_levels)) +
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
    for (i in seq_len(nrow(quantile_pairs))) {

      q_pair <- quantile_pairs[i, ]

      loop_txt <- paste0(
        "ggplot2::geom_ribbon(
          data    = crB_data,
          mapping = ggplot2::aes(x    = .data$dose,
                                 ymin = .data$'", as.character(q_pair[1]),"',
                                 ymax = .data$'", as.character(q_pair[2]),"'),
          fill    = acc_color,
          alpha   = 0.25)")

      plts <- plts + eval(parse(text = loop_txt))

    }
    rm(i)
    
    ## Bootstrapped Fit
    plts <- plts +
      ggplot2::geom_line(
        data    = crB_data,
        mapping = ggplot2::aes(.data$dose, .data$`0.5`),
        color   = acc_color)

  }

  plts <- plts +
    ## Posterior Credible Intervals
    {if (cr_intv) {

      ggplot2::geom_errorbar(
        data    = data.frame(x    = dose_levels,
                             ymin = post_quantiles[1, ],
                             ymax = post_quantiles[3, ]),
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
        medians     = post_quantiles[2, ]),
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