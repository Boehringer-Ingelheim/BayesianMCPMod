#' @title getModelFits
#' 
#' @description Fits dose-response curves for the specified dose-response models, based on the posterior distributions.
#' For the simplified fit, multivariate normal distributions will be approximated and reduced by one-dimensional normal distributions. 
#' For the default case, the Nelder-Mead algorithm is used. 
#' In detail, for both approaches the mean vector \eqn{\theta^{Y}} and the covariance \eqn{\Sigma} of the (mixture) posterior distributions and the corresponding posterior weights \eqn{\tilde{\omega}_{l}} for \eqn{l \in {1,...,L}} are used as basis
#' For the full fit a GLS estimator is used to minimize the following expression for the respective dose-response models \eqn{m}
#' \deqn{ \hat{\theta}_{m}=argmin_{\theta_{m}} \sum_{l=1}^{L} \tilde{\omega}_{l}(\theta_{l_{i}}^{Y}-f(dose_{i},\hat{\theta}_{m}))'\Sigma_{l}^{-1}(\theta_{l_{i}}^{Y}-f(dose_{i},\hat{\theta}_{m}))}
#' Therefore the function nloptr of the nloptr package is utilized. 
#' In the simplified case \eqn{L=1}, as the dimension of the posterior is reduced to 1 first.
#' The generalized AIC values are calculated via the formula
#' \deqn{gAIC_{m} = \sum_{l=1}^{L} \tilde{\omega}_{l} \sum_{i=0}^{K} \frac{1}{\Sigma_{l_{i,i}}} (\theta_{l_i}^Y - f(dose_{i},\hat{\theta}_{m}))^2 + 2p }
#' where \eqn{p} denotes the number of estimated parameters and \eqn{K} the number of active dose levels.
#' Here as well for the simplified case the formula reduces to one summand as \eqn{L=1}.
#' Corresponding gAIC based weights for model \eqn{M} are calculated as outlined in Schorning et al. (2016)
#' \deqn{
#' \Omega_I (M) = \frac{\exp(-0.5 gAIC_{M})}{\sum_{m=1}^{Q} \exp(-0.5 gAIC_{m})}
#' }
#' where \eqn{Q} denotes the number of models included in the averaging procedure.
#' @references Schorning K, Bornkamp B, Bretz F, Dette H. 2016. Model selection versus model averaging in dose finding studies. Stat Med; 35; 4021-4040.
#' @param models List of model names for which a fit will be performed.
#' @param dose_levels A vector containing the different dosage levels.
#' @param posterior A getPosterior object, containing the (multivariate) posterior distribution per dosage level. 
#' @param simple Boolean variable, defining whether simplified fit will be applied. Default FALSE.
#' @examples
#' posterior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 1), sigma = 2),
#'                        DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 3, s = 1.2), sigma = 2),
#'                        DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 4, s = 1.5), sigma = 2) ,  
#'                        DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 6, s = 1.2), sigma = 2) ,
#'                        DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 6.5, s = 1.1), sigma = 2))
#' models         <- c("emax", "exponential", "sigEmax", "linear")
#' dose_levels    <- c(0, 1, 2, 4, 8)
#' 
#' fit        <- getModelFits(models      = models,
#'                            posterior   = posterior_list,
#'                            dose_levels = dose_levels)
#' fit_simple <- getModelFits(models      = models,
#'                            posterior   = posterior_list,
#'                            dose_levels = dose_levels,
#'                            simple      = TRUE)
#'                            
#' @return An object of class modelFits. A list containing information about the fitted model coefficients, the prediction per dose group as well as maximum effect and generalized AIC (and corresponding weight) per model.
#' 
#' @export
getModelFits <- function (
    
  models,
  dose_levels,
  posterior,
  simple = FALSE
  
) {
  
  checkmate::check_list(models, any.missing = FALSE)
  checkmate::check_double(dose_levels, lower = 0, any.missing = FALSE, len = length(models))
  checkmate::check_class(posterior, "postList")
  checkmate::check_logical(simple)
  
  models      <- unique(gsub("\\d", "", models))
  
  getModelFit <- ifelse(simple, getModelFitSimple, getModelFitOpt)
  model_fits  <- lapply(models, getModelFit, dose_levels, posterior)
  
  model_fits  <- addModelWeights(model_fits)
  
  names(model_fits)             <- models
  attr(model_fits, "posterior") <- posterior
  class(model_fits)             <- "modelFits"
  
 
  return (model_fits)
  
}

## ignoring mixture posterior
getModelFitSimple <- function (
    
  model,
  dose_levels,
  posterior
  
) {
  
  fit <- DoseFinding::fitMod(
    dose  = dose_levels,
    resp  = summary.postList(posterior)[, 1],
    S     = diag(summary.postList(posterior)[, 2]^2),
    model = model,
    type  = "general",
    bnds  = DoseFinding::defBnds(mD = max(dose_levels))[[model]])
  
  pred_vals  <- stats::predict(fit, predType = "ls-means")
  max_effect <- max(pred_vals) - min(pred_vals)
  gAIC       <- DoseFinding::gAIC(fit)
  
  model_fit <- list(
    model       = model,
    coeffs      = fit$coefs,
    dose_levels = dose_levels,
    pred_values = pred_vals,
    max_effect  = max_effect,
    gAIC        = gAIC)
  
  return (model_fit)
  
}

## taking mixture posterior into account
getModelFitOpt <- function (
    
  model,
  dose_levels,
  posterior
  
) {
  
  getOptParams <- function (
    
    model,
    dose_levels,
    posterior
    
  ) {
    
    switch (model,
            "emax" = {
              lb     <- c(-Inf, -Inf, 0.001 * max(dose_levels))
              ub     <- c(Inf, Inf, 1.5 * max(dose_levels))
              expr_i <- quote(sum((post_combs$means[i, ] - (theta[1] + (theta[2] * dose_levels) / (theta[3] + dose_levels)))^2 / (post_combs$vars[i, ])))},
            "sigEmax" = {
              lb     <- c(-Inf, -Inf, 0.001 * max(dose_levels), 0.5)
              ub     <- c(Inf, Inf, 1.5 * max(dose_levels), 10)
              expr_i <- quote(sum((post_combs$means[i, ] - (theta[1] + (theta[2] * dose_levels^theta[4]) / (theta[3]^theta[4] + dose_levels^theta[4])))^2 / (post_combs$vars[i, ])))},
            "exponential" = {
              lb     <- c(-Inf, -Inf, 0.1 * max(dose_levels))
              ub     <- c(Inf, Inf, 2 * max(dose_levels))
              expr_i <- quote(sum((post_combs$means[i, ] - (theta[1] + theta[2] * (exp(dose_levels / theta[3]) - 1)))^2 / (post_combs$vars[i, ])))},
            "quadratic" = {
              lb     <- NULL
              ub     <- c(Inf, Inf, 0)
              expr_i <- quote(sum((post_combs$means[i, ] - (theta[1] + theta[2] * dose_levels + theta[3] * dose_levels^2))^2 / (post_combs$vars[i, ])))},
            "linear" = {
              lb     <- NULL
              ub     <- NULL
              expr_i <- quote(sum((post_combs$means[i, ] - (theta[1] + theta[2] * dose_levels))^2 / (post_combs$vars[i, ])))},
            "logistic" = {
              lb     <- c(-Inf, -Inf, 0.001 * max(dose_levels), 0.01 * max(dose_levels))
              ub     <- c(Inf, Inf, 1.5 * max(dose_levels), 0.5 * max(dose_levels))
              expr_i <- quote(sum((post_combs$means[i, ] - (theta[1] + theta[2] / (1 + exp((theta[3] - dose_levels) / theta[4]))))^2 / (post_combs$vars[i, ])))},
            {
              stop ("error")})
    
    simple_fit <- getModelFitSimple(
      model       = model,
      dose_levels = dose_levels,
      posterior   = posterior)
    
    param_list <- list(
      expr_i  = expr_i,
      opts    = list("algorithm" = "NLOPT_LN_NELDERMEAD", maxeval = 1e3),
      params  = list(
        x0 = simple_fit$coeffs,
        lb = lb,
        ub = ub),
      c_names = names(simple_fit$coeffs))
    
    return (param_list)
    
  }
  
  optFun <- function (
    
    theta,
    
    dose_levels,
    post_combs,
    expr_i
    
  ) {
    
    comps <- sapply(seq_along(post_combs$weights), function (i) eval(expr_i))
    
    return (sum(post_combs$weights * comps))
    
  }
  
  param_list <- getOptParams(model, dose_levels, posterior)
  post_combs <- getPostCombsI(posterior)
  
  fit <- nloptr::nloptr(
    eval_f = optFun,
    x0     = param_list$params$x0,
    lb     = param_list$params$lb,
    ub     = param_list$params$ub,
    opts   = param_list$opts,
    dose_levels = dose_levels,
    post_combs  = post_combs,
    expr_i      = param_list$expr_i)
  
  names(fit$solution) <- param_list$c_names
  
  model_fit <- list(
    model       = model,
    coeffs      = fit$solution,
    dose_levels = dose_levels)
  
  model_fit$pred_values <- predictModelFit(model_fit)
  model_fit$max_effect  <- max(model_fit$pred_values) - min(model_fit$pred_values)
  model_fit$gAIC        <- getGenAIC(model_fit, post_combs)
  
  return (model_fit)
  
}

predictModelFit <- function (
    
  model_fit,
  doses = NULL
  
) {
  
  if (is.null(doses)) {
    
    doses <- model_fit$dose_levels
    
  }
  
  pred_vals <- switch (
    model_fit$model,
    "emax" = {DoseFinding::emax(
      doses,
      model_fit$coeffs["e0"],
      model_fit$coeffs["eMax"],
      model_fit$coeffs["ed50"])},
    "sigEmax" = {DoseFinding::sigEmax(
      doses,
      model_fit$coeffs["e0"],
      model_fit$coeffs["eMax"],
      model_fit$coeffs["ed50"],
      model_fit$coeffs["h"])},
    "exponential" = {DoseFinding::exponential(
      doses,
      model_fit$coeffs["e0"],
      model_fit$coeffs["e1"],
      model_fit$coeffs["delta"])},
    "quadratic" = {DoseFinding::quadratic(
      doses,
      model_fit$coeffs["e0"],
      model_fit$coeffs["b1"],
      model_fit$coeffs["b2"])},
    "linear" = {DoseFinding::linear(
      doses,
      model_fit$coeffs["e0"],
      model_fit$coeffs["delta"])},
    "logistic" = {DoseFinding::logistic(
      doses,
      model_fit$coeffs["e0"],
      model_fit$coeffs["eMax"],
      model_fit$coeffs["ed50"],
      model_fit$coeffs["delta"])},
    {stop("error")})
  
  return (pred_vals)
  
}

getGenAIC <- function (
    
  model_fit,
  post_combs
  
) {
  
  expr_i   <- quote(sum(1 / post_combs$vars[i, ] * (post_combs$means[i, ] - model_fit$pred_values)^2) + 2 * length(model_fit$coeffs))
  comb_aic <- sapply(seq_along(post_combs$weights), function (i) eval(expr_i))
  
  g_aic    <- stats::weighted.mean(comb_aic, w = post_combs$weights)
  
  return (g_aic)
  
}

addModelWeights <- function (
    
  model_fits
  
) {
  
  g_aics        <- sapply(model_fits, function (model_fit) model_fit$gAIC)
  exp_values    <- exp(-0.5 * g_aics)
  model_weights <- exp_values / sum(exp_values)
  
  names(model_weights) <- NULL
  
  model_fits_out <- lapply(seq_along(model_fits), function (i) {
    
    c(model_fits[[i]], model_weight = model_weights[i])
    
  })
  
  return (model_fits_out)
  
}