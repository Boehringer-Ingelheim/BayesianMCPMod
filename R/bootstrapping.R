#' @title getBootstrapQuantiles
#'
#' @description A function for the calculation of credible intervals to assess the uncertainty for the model fit. 
#' Hereby the credible intervals are calculated as follows.
#' Samples from the posterior distribution are drawn (via the RBesT function rmix) and for every sample the simplified fitting step (see getModelFits function) and a prediction is performed. 
#' These fits are then used to identify the specified quantiles. 
#' This approach can be considered as the bayesian equivalent of the frequentist bootstrap approach described in O'Quigley, Iasonos, J., and B. Bornkamp. 2017. Handbook of Methods for Designing, Monitoring, and Analyzing Dose-Finding Trials. Boca Raton: CRC Press. 
#' Instead of drawing n bootstrap samples from the sampling distribution of the trial dose-response estimates, here the samples are directly taken from the posterior distribution.
#' 
#' 
#' @param model_fits an object of class modelFits, i.e. information about fitted models & corresponding model coefficients as well as the posterior distribution that was the basis for the model fitting 
#' @param quantiles a vector of quantiles that should be evaluated 
#' @param n_samples number of samples that should be drawn as basis for the 
#' @param doses a vector of doses for which a prediction should be performed
#' @param avg_fit boolean variable, defining whether an average fit (based on generalized AIC weights) should be performed in addition to the individual models. Default TRUE.
#'
#' @return  A data frame with entries doses, models, and quantiles 
#' @export
getBootstrapQuantiles <- function (

  model_fits,
  quantiles,
  n_samples = 1e3,
  doses     = NULL,
  avg_fit   = TRUE

) {
  
  mu_hat_samples <- sapply(attr(model_fits, "posterior"),
                           RBesT::rmix, n = n_samples)
  sd_hat         <- summary.postList(attr(model_fits, "posterior"))[, 2]
  
  dose_levels    <- model_fits[[1]]$dose_levels
  model_names    <- names(model_fits)
  quantile_probs <- sort(unique(quantiles))
  
  if (is.null(doses)) {
    
    doses <- seq(min(dose_levels), max(dose_levels), length.out = 100L)
    
  }
  
  preds <- apply(mu_hat_samples, 1, function (mu_hat) {
    
    preds_mu_hat <- sapply(model_names, function (model) {
      
      fit <- DoseFinding::fitMod(
        dose  = model_fits[[1]]$dose_levels,
        resp  = mu_hat,
        S     = diag(sd_hat^2),
        model = model,
        type  = "general",
        bnds  = DoseFinding::defBnds(
          mD = max(model_fits[[1]]$dose_levels))[[model]])
      
      preds <- stats::predict(fit, doseSeq = doses, predType = "ls-means")
      attr(preds, "gAIC") <- DoseFinding::gAIC(fit)
      
      return (preds)
      
    }, simplify = FALSE)
    
    preds_mu_mat <- do.call(rbind, preds_mu_hat)
    
    if (avg_fit) {
      
      avg_fit_indx <- which.min(sapply(preds_mu_hat, attr, "gAIC"))
      preds_mu_mat <- rbind(preds_mu_mat, avgFit = preds_mu_mat[avg_fit_indx, ])
      
    }
    
    return (preds_mu_mat)
    
  })
  
  if (avg_fit) {
    
    model_names <- c(model_names, "avgFit")
    
  }
  
  sort_indx <- order(rep(seq_along(model_names), length(doses)))
  quant_mat <- t(apply(X      = preds[sort_indx, ],
                       MARGIN = 1,
                       FUN    = stats::quantile,
                       probs  = quantile_probs))
  
  cr_bounds_data <- cbind(
    doses  = doses,
    models = rep(
      x    = factor(model_names,
                    levels = c("linear", "emax", "exponential",
                               "sigEmax", "logistic", "quadratic",
                               "avgFit")),
      each = length(doses)),
    as.data.frame(quant_mat))

  return (cr_bounds_data)

}
