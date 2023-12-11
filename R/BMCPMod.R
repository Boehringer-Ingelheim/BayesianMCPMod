#' @title assessDesign
#' 
#' @param n_patients tbd
#' @param mods tbd
#' @param prior_list tbd
#' @param n_sim tbd
#' @param alpha_crit_val tbd
#' @param simple tbd
#' 
#' @export
assessDesign <- function (
    
  n_patients,
  mods,
  prior_list,
  
  sd             = NULL,
  
  n_sim          = 1e3,
  alpha_crit_val = 0.05,
  simple         = TRUE
  
) {
  
  dose_levels <- attr(prior_list, "dose_levels")
  sd          <- ifelse(is.null(sd), attr(prior_list, "sd_tot"), sd)
  
  stopifnot(
    "sd length must coincide with number of dose levels" =
      length(sd) == length(dose_levels))
  
  data <- simulateData(
    n_patients  = n_patients,
    dose_levels = dose_levels,
    sd          = sd,
    mods        = mods,
    n_sim       = n_sim)
  
  model_names <- names(mods)
  
  eval_design <- lapply(model_names, function (model_name) {
    
    posterior_list <- getPosterior(
      data       = getModelData(data, model_name),
      prior_list = prior_list)
    
    crit_prob_adj <- getCritProb(
      mods           = mods,
      dose_levels    = dose_levels,
      dose_weights   = n_patients,
      alpha_crit_val = alpha_crit_val)
    
    contr_mat_prior <- getContr(
      mods           = mods,
      dose_levels    = dose_levels,
      dose_weights   = n_patients,
      prior_list     = prior_list)
    
    b_mcp_mod <- performBayesianMCPMod(
      posteriors_list = posterior_list,
      contr_mat       = contr_mat_prior,
      crit_prob_adj   = crit_prob_adj,
      simple          = simple)
    
  })
  
  names(eval_design) <- model_names
  
  return (eval_design)
  
}

#' @title getContr
#' 
#' @param mods tbd
#' @param dose_levels tbd
#' @param dose_weights tbd
#' @param prior_list tbd
#' @param se_new_trial tbd
#' @param sd_posterior tbd
#' 
#' @export
getContr <- function (
    
  mods,
  dose_levels,
  dose_weights = NULL,
  prior_list   = NULL,
  se_new_trial = NULL,
  sd_posterior = NULL
  
) {
  
  # frequentist & re-estimation
  if (!is.null(se_new_trial) & 
      is.null(dose_weights) & is.null(prior_list) & is.null(sd_posterior)) {
    
    w <- NULL
    S <- diag((se_new_trial)^2)
    
  # frequentist & no re-estimation
  } else if (!is.null(dose_weights) & 
             is.null(se_new_trial) & is.null(prior_list) & is.null(sd_posterior)) {
    
    w <- dose_weights
    S <- NULL
    
  # Bayesian & re-estimation
  } else if (!is.null(sd_posterior) & 
             is.null(se_new_trial) & is.null(prior_list) & is.null(dose_weights)) {
    
    w <- NULL
    S <- diag((sd_posterior)^2)
    
  # Bayesian & no re-estimation
  } else if (!is.null(dose_weights) & !is.null(prior_list) & 
             is.null(se_new_trial) & is.null(sd_posterior)) {
    
    w <- dose_weights +
      suppressMessages(round(unlist(lapply(prior_list, RBesT::ess))))
    S <- NULL
    
  } else {
    
    stop (paste("Provided combiations of 'se_new_trial',",
                 "'dose_weights', 'prior_list', 'sd_posterior' not allowed.",
                 "See ?getContr for allowed combinations."))
    
  }
  
  if (is.null(w)) {
    
    contr <- DoseFinding::optContr(
      models = mods,
      doses  = dose_levels,
      S      = S)
    
  } else {
    
    contr <- DoseFinding::optContr(
      models = mods,
      doses  = dose_levels,
      w      = w)
    
  }
  
  return (contr)
  
}

#' @title getCritProb
#' 
#' @param mods tbd
#' @param dose_levels tbd
#' @param dose_weights tbd
#' @param alpha_crit_val tbd
#' 
#' @export
getCritProb <- function (
    
  mods,
  dose_levels,
  dose_weights,
  alpha_crit_val
  
) {
  
  contr_mat <- DoseFinding::optContr(
    models = mods,
    doses  = dose_levels,
    w      = dose_weights)
  
  crit_prob <- pnorm(DoseFinding:::critVal(
    corMat      = contr_mat$corMat,
    alpha       = alpha_crit_val,
    df          = 0,
    alternative = "one.sided"))
  
  return (crit_prob)
  
}

#' @title performBayesianMCPMod
#' 
#' @param posteriors_list tbd
#' @param contr_mat tbd
#' @param crit_prob_adj tbd
#' @param simple tbd
#' 
#' @export
performBayesianMCPMod <- function (
    
  posteriors_list,
  contr_mat,
  crit_prob_adj,
  simple = FALSE
  
) {
  
  if (class(posteriors_list) == "postList") {
    
    posteriors_list <- list(posteriors_list)
    
  }
  
  b_mcp <- performBayesianMCP(
    posteriors_list = posteriors_list,
    contr_mat       = contr_mat,
    crit_prob_adj       = crit_prob_adj)
  
  model_shapes <- colnames(contr_mat$contMat)
  dose_levels  <- as.numeric(rownames(contr_mat$contMat))
  
  fits_list <- lapply(seq_along(posteriors_list), function (i) {
    
    if (b_mcp[i, 1]) {
      
      sign_models <- b_mcp[i, -c(1, 2)] > attr(b_mcp, "crit_prob_adj")
      
      model_fits  <- getModelFits(
        models      = model_shapes,
        dose_levels = dose_levels,
        posterior   = posteriors_list[[i]],
        simple      = simple)
      
      model_fits  <- addSignificance(model_fits, sign_models)
      
    } else {
      
      NULL
      
    }
    
  })
  
  bmcpmod        <- list(BayesianMCP = b_mcp, Mod = fits_list)
  class(bmcpmod) <- "BayesianMCPMod"
  
  return (bmcpmod)
  
}

addSignificance <- function (
    
  model_fits,
  sign_models

) {
  
  names(sign_models) <- NULL
  
  model_fits_out <- lapply(seq_along(model_fits), function (i) {
    
    c(model_fits[[i]], significant = sign_models[i])
    
  })
  
  attributes(model_fits_out) <- attributes(model_fits)
  
  return (model_fits_out)
  
}

#' @title performBayesianMCP
#' 
#' @param posteriors_list tbd
#' @param contr_mat tbd
#' @param crit_prob_adj tbd
#' 
#' @export
performBayesianMCP <- function(
    
  posteriors_list,
  contr_mat,
  crit_prob_adj
  
) {
  
  if (class(posteriors_list) == "postList") {
    
    posteriors_list <- list(posteriors_list)
    
  }
  
  b_mcp <- t(sapply(posteriors_list, BayesMCPi, contr_mat, crit_prob_adj))
  
  attr(b_mcp, "crit_prob_adj") <- crit_prob_adj
  class(b_mcp) <- "BayesianMCP"
  
  return (b_mcp)
  
}

BayesMCPi <- function (
    
  posterior_i,
  contr_mat,
  crit_prob_adj
  
) {
  
  getPostProb <- function (
    
    contr_j,     # j: dose level
    post_combs_i # i: simulation outcome
    
  ) {
    
    ## Test statistic = sum over all components of
    ## posterior weight * normal probability distribution of
    ## critical values for doses * estimated mean / sqrt(product of critical values for doses)
    
    ## Calculation for each component of the posterior
    contr_theta   <- apply(post_combs_i$means, 1, `%*%`, contr_j)
    contr_var     <- apply(post_combs_i$vars, 1, `%*%`, contr_j^2)
    contr_weights <- post_combs_i$weights
    
    ## P(c_m * theta > 0 | Y = y) for a shape m (and dose j)
    post_probs <- sum(contr_weights * stats::pnorm(contr_theta / sqrt(contr_var)))
    
    return (post_probs)
    
  }
  
  post_combs_i <- getPostCombsI(posterior_i)
  post_probs   <- apply(contr_mat$contMat, 2, getPostProb, post_combs_i)
  
  res <- c(sign          = ifelse(max(post_probs) > crit_prob_adj, 1, 0),
           crit_prob_adj = crit_prob_adj,
           max_post_prob = max(post_probs),
           post_probs    = post_probs)
  
  return (res)
  
}
