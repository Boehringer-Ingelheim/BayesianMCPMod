#' @title getBootstrapBands
#'
#' @param model_fits tbd
#' @param n_samples tbd
#' @param alpha tbd
#' @param avg_fit tbd
#' @param dose_seq tbd
#'
#' @return tbd
#' @export
getBootstrapBands <- function (

  model_fits,
  n_samples = 1e3,
  alpha     = c(0.05, 0.5),
  avg_fit   = TRUE,
  dose_seq  = NULL

) {
  
  mu_hat_samples <- sapply(attr(model_fits, "posterior"),
                           RBesT::rmix, n = n_samples)
  sd_hat         <- summary.postList(attr(model_fits, "posterior"))[, 2]
  
  dose_levels    <- model_fits[[1]]$dose_levels
  model_names    <- names(model_fits)
  quantile_probs <- c(0.5, sort(unique(c(alpha / 2, 1 - alpha / 2))))
  
  if (is.null(dose_seq)) {
    
    dose_seq <- seq(min(dose_levels), max(dose_levels), length.out = 100L)
    
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
      
      preds <- stats::predict(fit, doseSeq = dose_seq, predType = "ls-means")
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
  
  sort_indx <- order(rep(seq_along(model_names), length(dose_seq)))
  quant_mat <- t(apply(X      = preds[sort_indx, ],
                       MARGIN = 1,
                       FUN    = stats::quantile,
                       probs  = quantile_probs))
  
  cr_bounds_data <- cbind(
    dose_seqs = dose_seq,
    models    = rep(
      x    = factor(model_names,
                    levels = c("linear", "emax", "exponential",
                               "sigEmax", "logistic", "quadratic",
                               "avgFit")),
      each = length(dose_seq)),
    as.data.frame(quant_mat))

  return (cr_bounds_data)

}
