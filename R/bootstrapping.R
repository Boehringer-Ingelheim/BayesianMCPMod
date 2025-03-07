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
#' @return  A data frame with columns for model, dose, and bootstrapped samples
#'
#' @export
getBootstrapQuantiles <- function (

  model_fits,
  quantiles,
  n_samples = 1e3,
  doses     = NULL

) {
  
  ## R CMD --as-cran appeasement
  models <- q_probs <- NULL

  bs_samples <- getBootstrapSamples(
    model_fits = model_fits,
    n_samples  = n_samples,
    doses      = doses)

  quantile_probs <- sort(unique(quantiles))
    
  bs_quantiles <- bs_samples |>
    dplyr::group_by(models, doses) |>
    dplyr::summarize(
      quantiles = list(stats::quantile(sample, probs = quantile_probs)),
      .groups = 'drop') |>
    tidyr::unnest_wider(quantiles, names_sep = "_") |>
    dplyr::rename_with(~ paste0(quantile_probs, "%"), dplyr::starts_with("quantiles_")) |>
    tidyr::pivot_longer(cols      = tidyr::ends_with("%"),
                        names_to  = "q_probs",
                        values_to = "q_values") |>
    dplyr::mutate(q_probs = as.numeric(gsub("%", "", q_probs)))

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
#' @return  A data frame with entries model, dose, and sample
#' 
#' @export
getBootstrapSamples <- function (
    
  model_fits,
  n_samples = 1e3,
  doses     = NULL
  
) {
  
  ## R CMD --as-cran appeasement
  name <- NULL
  
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
  preds <- t(future.apply::future_apply(mu_hat_samples, 1, function (mu_hat) {
    
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
      preds_mu_mat <- rbind(preds_mu_mat, avgFit = preds_mu_mat[avg_fit_indx, ])
      
    }
    
    # predictions[models, doses]
    return (preds_mu_mat)
    
  }))
  
  if (avg_fit) {
    
    model_names <- c(model_names, "avgFit")
    
  }
  
  bs_samples           <- as.data.frame(preds)
  colnames(bs_samples) <- paste0(model_names, "_", rep(doses, each = length(model_names)))
  
  bs_samples_long <- bs_samples |>
    tidyr::pivot_longer(cols      = tidyr::everything(),
                        values_to = "sample") |>
    tidyr::separate(col  = name,
                    into = c("models", "doses"),
                    sep  = "_") |>
    dplyr::mutate(doses = as.numeric(doses))
  
  return (bs_samples_long)
  
}
