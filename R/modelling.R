#' @title getModelFits
#' 
#' @description Fits dose-response curves for the specified dose-repsonse models, based on the posterior distributions.
#' For the simplified fit, multivariate normal distributions will be approximated and reduced by one-dimensional normal distributions. 
#' For the default case, the Nelder-Mead algorithm is used. Will be further updated and links to publication as well as references will be added.
#' 
#' @param models list of model names for which a fit will be performed.
#' @param dose_levels a vector containing the different dosage levels.
#' @param posterior a getPosterior object, containing the (multivariate) posterior distribution per dosage level. 
#' @param simple boolean variable, defining whether simplified fit will be applied. Default FALSE.
#' 
#' @return model_fits returns a list, containing information about the fitted model coefficients, the prediction per dose group as well as maximum effect and generalized AIC per model.
#' 
#' @export
getModelFits <- function (
    
  models,
  dose_levels,
  posterior,
  simple = FALSE
  
) {
  
<<<<<<< HEAD
  checkmate::check_list(models, any.missing = FALSE)
  checkmate::check_double(dose_levels, lower = 0, any.missing = FALSE, len = length(models))
  checkmate::check_class(posterior, "postList")
  checkmate::check_logical(simple)
=======
  models      <- unique(gsub("\\d", "", models))
>>>>>>> 25f8f28541ce5be0c55fecee4cdd47a6c8603237
  
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
