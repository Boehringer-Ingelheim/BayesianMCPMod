#' @title getModelFits
#'
#' @description Fits dose-response curves for the specified dose-response models, based on the posterior distributions.
#' For the simplified fit, multivariate normal distributions will be approximated and reduced by one-dimensional normal distributions.
#' For the default case, the Nelder-Mead algorithm is used.
#' In detail, for both approaches the mean vector \eqn{\theta^{Y}} and the covariance \eqn{\Sigma} of the (mixture) posterior distributions and the corresponding posterior weights \eqn{\tilde{\omega}_{l}} for \eqn{l \in {1,...,L}} are used as basis
#' For the full fit a GLS estimator is used to minimize the following expression for the respective dose-response models \eqn{m}
#' \deqn{ \hat{\theta}_{m}=\text{arg min}_{\theta_{m}} \sum_{l=1}^{L} \tilde{\omega}_{l}(\theta_{l_{i}}^{Y}-f(dose_{i},\hat{\theta}_{m}))'\Sigma_{l}^{-1}(\theta_{l_{i}}^{Y}-f(dose_{i},\hat{\theta}_{m}))}
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
#' @param models List (or vector) of model names for which a fit will be performed.
#' @param dose_levels A vector containing the different dosage levels.
#' @param posterior A getPosterior object, containing the (multivariate) posterior distribution per dosage level.
#' @param avg_fit Boolean variable, defining whether an average fit (based on generalized AIC weights) should be performed in addition to the individual models. Default TRUE.
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
#'                            
#' fit
#'                            
#' fit_simple <- getModelFits(models      = models,
#'                            posterior   = posterior_list,
#'                            dose_levels = dose_levels,
#'                            simple      = TRUE)
#'                            
#' fit_simple
#'
#' @return An object of class modelFits. A list containing information about the fitted model coefficients, the prediction per dose group as well as maximum effect and generalized AIC (and corresponding weight) per model.
#'
#' @export
getModelFits <- function (

  models,
  dose_levels,
  posterior,
  avg_fit = TRUE,
  simple  = FALSE

) {
  
  if (inherits(models, "character")) {
    models <- stats::setNames(as.list(models), models)
  }
  
  checkmate::check_list(models, any.missing = FALSE)
  checkmate::check_double(dose_levels, lower = 0, any.missing = FALSE, len = length(models))
  checkmate::check_class(posterior, "postList")
  checkmate::check_logical(avg_fit)
  checkmate::check_logical(simple)

  model_names <- unique(gsub("\\d", "", names(models)))

  getModelFit <- ifelse(simple, getModelFitSimple, getModelFitOpt)

  model_fits  <- lapply(model_names, getModelFit, dose_levels, posterior, list("scal" = attr(models, "scal")))
  model_fits  <- addModelWeights(model_fits)
  
  names(model_fits) <- model_names
  
  if (avg_fit) {
    
    model_fits <- addAvgFit(model_fits)
    
  }

  attr(model_fits, "posterior") <- posterior
  class(model_fits)             <- "modelFits"

  return (model_fits)

}

addAvgFit <- function (
    
  model_fits,
  doses = NULL
  
) {
  
  pred_vals  <- predictAvgFit(model_fits = model_fits, doses = doses)
  
  avg_fit  <- list(model        = "avgFit",
                   coeff        = NA,
                   dose_levels  = model_fits[[1]]$dose_levels,
                   pred_values  = pred_vals,
                   max_effect   = max(pred_vals) - min(pred_vals),
                   gAIC         = NA,
                   model_weight = NA)
  
  model_fits <- c(model_fits, avgFit = list(avg_fit))
  
  return (model_fits)
  
}

predictAvgFit <- function (
    
  model_fits,
  doses = NULL
    
) {
  
  model_fits_sub <- model_fits[names(model_fits) != "avgFit"]
  
  if (is.null(doses)) {
    
    dose_levels <- model_fits_sub[[1]]$dose_levels
    
    stopifnot(all(sapply(model_fits_sub,
                         function (x) identical(x$dose_levels, dose_levels))))
    
    mod_preds <- sapply(model_fits_sub, function (x) x$pred_values)
    
  } else {
    
    mod_preds <- do.call(cbind, predict.modelFits(model_fits_sub, doses = doses))
    
  }
  
  mod_weights <- sapply(model_fits_sub, function (x) x$model_weight)
  
  pred_vals   <- as.vector(mod_preds %*% mod_weights)
  
  return (pred_vals)
  
}

## ignoring mixture posterior
getModelFitSimple <- function (

  model,
  dose_levels,
  posterior,
  add_args = NULL

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
  posterior,
  add_args = NULL

) {

  getOptParams <- function (

    model,
    dose_levels,
    posterior,
    add_args = NULL

  ) {
    
    switch (
      model,
      "emax" = {
        lb     <- c(-Inf, -Inf, 0.001 * max(dose_levels))
        ub     <- c(Inf, Inf, 1.5 * max(dose_levels))
        expr_i <- quote(sum((post_combs$means[i, ] - (theta[1] + (theta[2] * dose_levels) / (theta[3] + dose_levels)))^2 / (post_combs$vars[i, ])))
      },
      "sigEmax" = {
        lb     <- c(-Inf, -Inf, 0.001 * max(dose_levels), 0.5)
        ub     <- c(Inf, Inf, 1.5 * max(dose_levels), 10)
        expr_i <- quote(sum((post_combs$means[i, ] - (theta[1] + (theta[2] * dose_levels^theta[4]) / (theta[3]^theta[4] + dose_levels^theta[4])))^2 / (post_combs$vars[i, ])))
      },
      "exponential" = {
        lb     <- c(-Inf, -Inf, 0.1 * max(dose_levels))
        ub     <- c(Inf, Inf, 2 * max(dose_levels))
        expr_i <- quote(sum((post_combs$means[i, ] - (theta[1] + theta[2] * (exp(dose_levels / theta[3]) - 1)))^2 / (post_combs$vars[i, ])))
      },
      "quadratic" = {
        lb     <- NULL
        ub     <- c(Inf, Inf, Inf)
        expr_i <- quote(sum((post_combs$means[i, ] - (theta[1] + theta[2] * dose_levels + theta[3] * dose_levels^2))^2 / (post_combs$vars[i, ])))
      },
      "linear" = {
        lb     <- NULL
        ub     <- NULL
        expr_i <- quote(sum((post_combs$means[i, ] - (theta[1] + theta[2] * dose_levels))^2 / (post_combs$vars[i, ])))
      },
      "logistic" = {
        lb     <- c(-Inf, -Inf, 0.001 * max(dose_levels), 0.01 * max(dose_levels))
        ub     <- c(Inf, Inf, 1.5 * max(dose_levels), 0.5 * max(dose_levels))
        expr_i <- quote(sum((post_combs$means[i, ] - (theta[1] + theta[2] / (1 + exp((theta[3] - dose_levels) / theta[4]))))^2 / (post_combs$vars[i, ])))
      },
      "betaMod" = {
        lb     <- c(-Inf, -Inf, 0.05, 0.05)
        ub     <- c(Inf, Inf, 4, 4)
        scal   <- ifelse(is.null(add_args), 1.2 * max(dose_levels), add_args[["scal"]]) #for betaMod shape
        expr_i <- substitute(sum(((post_combs$means[i, ] - (theta[1] + theta[2] * (((theta[3] + theta[4])^(theta[3] + theta[4])) / ((theta[3] ^ theta[3]) * (theta[4]^theta[4]))) * (dose_levels / scal)^(theta[3]) * ((1 - dose_levels / scal)^(theta[4]))))^2 / (post_combs$vars[i, ]))),
                             list(scal = scal))
      },
      stop(paste("error:", model, "shape not (yet) implemented"))
    )

    simple_fit <- getModelFitSimple(
      model       = model,
      dose_levels = dose_levels,
      posterior   = posterior)

    param_list <- list(
      expr_i  = expr_i,
      opts    = list("algorithm" = "NLOPT_LN_NELDERMEAD", maxeval = 1e3), #,stopval=0
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
    expr_i      = param_list$expr_i
  )

  names(fit$solution) <- param_list$c_names

  model_fit <- list(
    model       = model,
    coeffs      = fit$solution,
    dose_levels = dose_levels)

  model_fit$pred_values <- predictModelFit(model_fit = model_fit, doses = model_fit$dose_levels, add_args = add_args)
  model_fit$max_effect  <- max(model_fit$pred_values) - min(model_fit$pred_values)
  model_fit$gAIC        <- getGenAIC(model_fit, post_combs)

  return (model_fit)

}



predictModelFit <- function (

  model_fit,
  doses    = NULL,
  add_args = NULL

) {
  
  if (is.null(doses)) {
    doses <- model_fit$dose_levels
  }

  ## cannot predict average values as model_fits required instead of model_fit
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
    "betaMod" = {DoseFinding::betaMod(
      doses,
      model_fit$coeffs["e0"],
      model_fit$coeffs["eMax"],
      model_fit$coeffs["delta1"],
      model_fit$coeffs["delta2"],
      ifelse(is.null(add_args) | is.null(add_args[["scal"]]), 1.2 * max(doses), add_args[["scal"]]))},
    "avgFit" = {stop("error: use predictAvgFit for the avgFit")},
    {stop(paste("error:", model_fit$model, "shape not yet implemented"))})

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

#' @title getMED
#'
#' @description This function provides information on the minimally efficacious dose (MED).
#' The MED evaluation can either be based on the fitted model shapes (model_fits) or on bootstrapped quantiles (bs_quantiles). 
#' @details
#' The function assumes that the 1st dose group is the control dose group.
#' 
#' The bootstrap approach allows for an MED based on decision rules of the form
#' \deqn{\widehat{\text{MED}} = \text{arg min}_{d\in\{d_1, \dots, d_k\}} \left\{ \text{Pr}\left(f(d, \hat\theta) - f(d_1, \hat\theta) > \Delta\right) > \gamma \right\} .}
#' The model-shape approach takes the point estimate of the model into account.
#' @param delta A numeric value for the threshold Delta.
#' @param evidence_level A numeric value between 0 and 1 for the evidence level gamma. Used for the bs_quantiles-based evaluation and not used for the model_fits-based evaluation. Default 0.5.
#' @param dose_levels A vector of numerics containing the different dosage levels. Default NULL.
#' @param model_fits An object of class modelFits as created with getModelFits(). Default NULL.
#' @param bs_quantiles A dataframe created with getBootstrapQuantiles(). Default NULL.
#' @return  A matrix with rows for MED reached, MED, and MED index in the vector of dose levels and columns for the dose-response shapes.
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
#' # MED based on the model_fit:
#' getMED(delta = 5, model_fits = model_fits)
#'                                
#' # MED based on bootstrapped quantiles
#' bs_quantiles <- getBootstrapQuantiles(model_fits = model_fits,
#'                                       quantiles  = c(0.025, 0.2, 0.5),
#'                                       n_samples  = 100) # speeding up example run time
#'                                       
#' getMED(delta          = 5,
#'        evidence_level = 0.8,
#'        bs_quantiles   = bs_quantiles)
#'
#' @export
getMED <- function (
    
  delta,
  evidence_level = 0.5,
  dose_levels    = NULL,
  model_fits     = NULL,
  bs_quantiles   = NULL
  
) {
  
  stopifnot("Either model_fits or bs_quantiles must be not NULL, but not both" =
              is.null(model_fits) & !is.null(bs_quantiles) |
              !is.null(model_fits) & is.null(bs_quantiles))
  
  checkmate::check_double(delta)
  checkmate::check_double(evidence_level, lower = 0, upper = 1)
  checkmate::check_double(dose_levels, lower = 0, any.missing = FALSE)
  
  if (!is.null(model_fits)) {
    
    checkmate::check_class(model_fits, "modelFits")
    
    if (is.null(dose_levels)) {
      
      dose_levels <- model_fits[[1]]$dose_levels
      
    }
    
    preds <- do.call(rbind, predict.modelFits(model_fits, doses = dose_levels))
    
  } else {
    
    if (is.null(dose_levels)) {
      
      dose_levels <- unique(bs_quantiles$doses)
      
    }
    
    # R CMD Check Appeasement
    q_probs <- doses <- q_values <- models <- NULL
    
    stopifnot("corresponding quantile (i.e. 1 - evidence_level) not in bootstrapped quantiles matrix" = 
                evidence_level %in% (1 - bs_quantiles$q_probs))
    
    stopifnot("dose_levels not in bootstrapped quantiles matrix" = 
                all(dose_levels %in% bs_quantiles$doses))
    
    preds <- bs_quantiles |>
      dplyr::filter((1 - q_probs) %in% evidence_level) |>
      dplyr::filter(doses %in% dose_levels) |>
      tidyr::pivot_wider(names_from = doses, values_from = q_values)
    
    model_names <- preds$models
    
    preds <- preds |>
      dplyr::select(-models, -q_probs) |>
      as.matrix()
    
    rownames(preds) <- model_names
    
  }
  
  med_info <- apply(preds, 1, function (model_preds) {
    
    # assumes 1st DG to be control DG
    med_indx_m <- abs(model_preds[2:length(model_preds)] - model_preds[1]) > delta
    
    if (!any(med_indx_m)) {
      
      med_reached <- 0
      med         <- NA_real_
      
    } else {
      
      med_reached <- 1
      med_indx    <- min(which(med_indx_m)) + 1
      med         <- dose_levels[med_indx]
      
    }
    
    return (c(med_reached = med_reached,
              med         = med))
    
  })
  
  return (med_info)
  
}

