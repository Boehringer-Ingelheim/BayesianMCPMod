#' @title getBootstrapQuantiles
#'
#' @description A function for the calculation of bootstrapped model predictions.
#' Samples from the posterior distribution are drawn (via the RBesT function rmix()) and for every sample the simplified fitting step (see getModelFits() function) and a prediction is performed.
#' These fits are then used to identify the specified quantiles.
#' This approach can be considered as the Bayesian equivalent of the frequentist bootstrap approach described in O'Quigley et al. (2017).
#' Instead of drawing n bootstrap samples from the sampling distribution of the trial dose-response estimates, here the samples are directly taken from the posterior distribution.
#' @references O'Quigley J, Iasonos A, Bornkamp B. 2017. Handbook of Methods for Designing, Monitoring, and Analyzing Dose-Finding Trials (1st ed.). Chapman and Hall/CRC. doi:10.1201/9781315151984
#'
#' @param model_fits An object of class modelFits, i.e. information about fitted models & corresponding model coefficients as well as the posterior distribution that was the basis for the model fitting
#' @param quantiles A vector of quantiles that should be evaluated
#' @param n_samples Number of samples that should be drawn as basis for the bootstrapped quantiles
#' @param doses A vector of doses for which a prediction should be performed. If NULL, the dose levels of the model_fits will be used. Default NULL.
#' @param probability_scale A boolean variable to specify if the trial has a continuous or a binary outcome. Setting to TRUE will transform predictions from the logit scale to the probability scale, which can be desirable for a binary outcome. Default `attr(model_fits, "probability_scale")`.
#'
#' @examples
#' posterior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 1), sigma = 2),
#'                        DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 3, s = 1.2), sigma = 2),
#'                        DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 4, s = 1.5), sigma = 2) ,
#'                        DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 6, s = 1.2), sigma = 2) ,
#'                        DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 6.5, s = 1.1), sigma = 2))
#' models         <- c("exponential", "linear")
#' dose_levels    <- c(0, 1, 2, 4, 8)
#' model_fits     <- getModelFits(models      = models,
#'                                posterior   = posterior_list,
#'                                dose_levels = dose_levels,
#'                                simple      = TRUE)
#'
#' bs_quantiles <- getBootstrapQuantiles(model_fits = model_fits,
#'                                       quantiles  = c(0.025, 0.5, 0.8, 0.975),
#'                                       n_samples  = 10, # speeding up example run time
#'                                       doses      = c(0, 6, 8))
#'                       
#' bs_quantiles
#' @return  A tibble with columns for model, dose, and bootstrapped samples
#'
#' @export
getBootstrapQuantiles <- function (

  model_fits,
  quantiles,
  n_samples         = 1e3,
  doses             = NULL,
  probability_scale = attr(model_fits, "probability_scale")

) {
  
  checkmate::assert_flag(probability_scale, null.ok = TRUE)
  if (is.null(probability_scale)) probability_scale <- FALSE
  
  ## R CMD --as-cran appeasement
  model <- dose <- sample_q <- sample_diff_q <- q_prob <- NULL

  bs_samples <- getBootstrapSamples(
    model_fits        = model_fits,
    n_samples         = n_samples,
    doses             = doses,
    probability_scale = probability_scale)

  quantile_probs <- sort(unique(quantiles))
    
  bs_quantiles <- bs_samples |>
    dplyr::group_by(model, dose) |>
    dplyr::summarize(
      sample_q      = list(stats::quantile(abs, probs = quantile_probs)),
      sample_diff_q = list(stats::quantile(diff, probs = quantile_probs)),
      .groups = "drop") |>
    tidyr::unnest_wider(sample_q, names_sep = "_") |>
    tidyr::unnest_wider(sample_diff_q, names_sep = "_") |>
    dplyr::rename_with(~ paste0("sample_", quantile_probs, "%"), dplyr::starts_with("sample_q_")) |>
    dplyr::rename_with(~ paste0("sample_diff_", quantile_probs, "%"), dplyr::starts_with("sample_diff_q_")) |>
    tidyr::pivot_longer(
      cols          = dplyr::matches("^(sample|sample_diff)_\\d+(\\.\\d+)?%$"),
      names_to      = c("sample_type", "q_prob"),
      names_pattern = "^(sample(?:_diff)?)_(.+)%$",
      values_to     = "q_val"
    ) |>
    dplyr::mutate(q_prob      = as.numeric(gsub("%", "", q_prob)),
                  sample_type = dplyr::case_when(
                    sample_type == "sample" ~ "abs",
                    sample_type == "sample_diff" ~ "diff"))
  
  attr(bs_quantiles, "direction") <- attr(model_fits, "direction")
  attr(bs_quantiles, "probability_scale") <- probability_scale
  
  return (bs_quantiles)

}

#' @title getBootstrapSamples
#'
#' @description A function to return bootstrap samples from the fitted dose-response models.
#' Samples from the posterior distribution are drawn (via the RBesT function rmix()) and for every sample the simplified fitting step (see getModelFits() function) and a prediction is performed. 
#' These samples are returned by this function.   
#' This approach can be considered as the Bayesian equivalent of the frequentist bootstrap approach described in O'Quigley et al. (2017).
#' Instead of drawing n bootstrap samples from the sampling distribution of the trial dose-response estimates, here the samples are directly taken from the posterior distribution.
#' @references O'Quigley J, Iasonos A, Bornkamp B. 2017. Handbook of Methods for Designing, Monitoring, and Analyzing Dose-Finding Trials (1st ed.). Chapman and Hall/CRC. doi:10.1201/9781315151984
#' @param model_fits An object of class modelFits, i.e. information about fitted models & corresponding model coefficients as well as the posterior distribution that was the basis for the model fitting 
#' @param n_samples Number of samples that should be drawn
#' @param doses A vector of doses for which a prediction should be performed
#' @param probability_scale A boolean variable to specify if the trial has a continuous or a binary outcome. Setting to TRUE will transform predictions from the logit scale to the probability scale, which can be desirable for a binary outcome. Default `attr(model_fits, "probability_scale")`.
#'
#' @examples
#' posterior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 1), sigma = 2),
#'                        DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 3, s = 1.2), sigma = 2),
#'                        DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 4, s = 1.5), sigma = 2) ,  
#'                        DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 6, s = 1.2), sigma = 2) ,
#'                        DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 6.5, s = 1.1), sigma = 2))
#' models         <- c("exponential", "linear")
#' dose_levels    <- c(0, 1, 2, 4, 8)
#' model_fits     <- getModelFits(models      = models,
#'                                posterior   = posterior_list,
#'                                dose_levels = dose_levels,
#'                                simple      = TRUE)
#'                            
#' bs_samples <- getBootstrapSamples(model_fits = model_fits,
#'                                   n_samples  = 10, # speeding up example run time
#'                                   doses      = c(0, 6, 8))
#'                       
#' bs_samples
#' @return  A tibble with columns for sample_id, model, dose, sample, and sample_diff
#' 
#' @export
getBootstrapSamples <- function (
    
  model_fits,
  n_samples         = 1e3,
  doses             = NULL,
  probability_scale = attr(model_fits, "probability_scale")
  
) {
  
  checkmate::assert_flag(probability_scale, null.ok = TRUE)
  if (is.null(probability_scale)) probability_scale <- FALSE
  
  ## R CMD --as-cran appeasement
  name <- dose <- sample_type <- sample_id <- NULL
  
  ## Make parallel processing optional
  if (requireNamespace("future.apply", quietly = TRUE)) {
    optPar_apply <- future.apply::future_apply
  } else {
    optPar_apply <- apply
  }
  
  mu_hat_samples <- sapply(attr(model_fits, "posterior"),
                           RBesT::rmix, n = n_samples)
  sd_hat         <- summary.postList(attr(model_fits, "posterior"))[, 2]
  
  dose_levels    <- model_fits[[1]]$dose_levels
  model_names    <- names(model_fits)
  
  if (is.null(doses)) {
    
    doses <- dose_levels
    
  }
  
  avg_fit <- "avgFit" %in% model_names
  
  if (avg_fit) {
    
    model_names <- setdiff(model_names, "avgFit")
    
  }
  
  # predictions[samples, (dose1 model1, dose1 model2, ..., dose2 model1, ...)]
  preds <- t(optPar_apply(mu_hat_samples, 1, function (mu_hat) {
    
    preds_mu_hat <- sapply(model_names, function (model) {
      
      fit <- DoseFinding::fitMod(
        dose  = dose_levels,
        resp  = mu_hat,
        S     = diag(sd_hat^2),
        model = model,
        type  = "general",
        bnds  = DoseFinding::defBnds(mD = max(dose_levels))[[model]])
      
      preds <- stats::predict(fit, doseSeq = doses, predType = "ls-means")
      attr(preds, "gAIC") <- DoseFinding::gAIC(fit)
      
      return (preds)
      
    }, simplify = FALSE)
    
    preds_mu_mat <- do.call(rbind, preds_mu_hat)
    
    if (avg_fit) {
      
      avg_fit_indx <- which.min(sapply(preds_mu_hat, attr, "gAIC"))
      preds_mu_mat <- rbind(avgFit = preds_mu_mat[avg_fit_indx, ], preds_mu_mat)
      
    }
    
    if (probability_scale) {

      preds_mu_mat <- RBesT::inv_logit(preds_mu_mat)

    }
    
    preds_mu_mat_adj <- preds_mu_mat - preds_mu_mat[, 1]
    
    # predictions[models, c(doses, adj_doses)]
    return (c(preds_mu_mat, preds_mu_mat_adj))
    
  }))
  
  if (avg_fit) {
    
    model_names <- c("avgFit", model_names)
    
  }
  
  bs_samples           <- as.data.frame(preds)
  colnames(bs_samples) <- c(paste0(model_names, "_",
                                   rep(doses, each = length(model_names))),
                            paste0(model_names, "_diff_",
                                   rep(doses, each = length(model_names))))
  
  bs_samples <- cbind(sample_id = seq_len(nrow(bs_samples)), bs_samples)
  
  bs_samples_long <- bs_samples |>
    tidyr::pivot_longer(
      cols          = -sample_id,
      names_to      = c("model", "diff", "dose"),
      names_pattern = "^(avgFit|linear|exponential|logistic|emax|sigEmax|quadratic|betaMod)(_diff)?_([0-9]+(?:\\.[0-9]+)?)$",
      values_to     = "sample") |>
    dplyr::mutate(
      diff = dplyr::case_when(
        diff == "" ~ "abs",
        diff == "_diff" ~ "diff"),
      dose = as.numeric(dose)) |>
    dplyr::rename(sample_type = diff) |>
    tidyr::pivot_wider(
      names_from  = sample_type,
      values_from = sample)
  
  attr(bs_samples_long, "probability_scale") <- probability_scale
  
  return (bs_samples_long)
  
}
