#' @title assessDesign
#'.
#' @description This function performs simulation based trial design evaluations for a set of specified dose-response models
#'
#' @param n_patients Vector specifying the planned number of patients per dose group
#' @param mods An object of class "Mods" as specified in the Dosefinding package.
#' @param prior_list a prior_list object specifying the utilized prior for the different dose groups 
#' @param sd tbd
#' @param n_sim number of simulations to be performed
#' @param alpha_crit_val critical value to be used for the testing (on the probability scale)
#' @param simple boolean variable, defining whether simplified fit will be applied. Passed to the getModelFits function. Default FALSE.
#' @param reestimate tbd Default FALSE
#' @param contr tbd Default NULL
#' @param dr_means tbd Default NULL
#' 
#' @export
assessDesign <- function (
    
  n_patients,
  mods,
  prior_list,
  
  sd,
  
  n_sim          = 1e3,
  alpha_crit_val = 0.05,
  simple         = TRUE,
  reestimate     = FALSE,
  
  contr          = NULL,
  dr_means       = NULL
  
) {
  
  dose_levels <- attr(mods, "doses")
  
  data <- simulateData(
    n_patients  = n_patients,
    dose_levels = dose_levels,
    sd          = sd,
    mods        = mods,
    n_sim       = n_sim, 
    dr_means    = dr_means)
  
  model_names <- colnames(data)[-c(1:3)]
  
  crit_prob_adj <- getCritProb(
    mods           = mods,
    dose_levels    = dose_levels,
    dose_weights   = n_patients,
    alpha_crit_val = alpha_crit_val)
  
  if (!reestimate & is.null(contr)) {
    
    contr <- getContr(
      mods           = mods,
      dose_levels    = dose_levels,
      dose_weights   = n_patients,
      prior_list     = prior_list)
    
  }
  
  eval_design <- lapply(model_names, function (model_name) {
    
    posterior_list <- getPosterior(
      data       = getModelData(data, model_name),
      prior_list = prior_list)
    
    if (reestimate & is.null(contr)) {
      
      post_sds <- sapply(posterior_list, function (post) summary(post)[, 2])
      
      contr <- apply(post_sds, 2, function (post_sd) getContr(
          mods           = mods,
          dose_levels    = dose_levels,
          sd_posterior   = post_sd))
      
    } 
    
    b_mcp_mod <- performBayesianMCPMod(
      posterior_list = posterior_list,
      contr          = contr,
      crit_prob_adj  = crit_prob_adj,
      simple         = simple)
    
  })
  
  names(eval_design) <- model_names
  
  attr(eval_design, "placEff")    <- attr(mods, "placEff")
  attr(eval_design, "maxEff")     <- attr(mods, "maxEff")
  attr(eval_design, "sampleSize") <- n_patients
  attr(eval_design, "priorESS")   <- getESS(prior_list)
  
  return (eval_design)
  
}

#' @title getContr
#' 
#' @description This function calculates contrast vectors that are optimal for detecting certain alternatives. More information and link to publication will be added.
#' 
#' @param mods An object of class "Mods" as specified in the Dosefinding package.
#' @param dose_levels vector containing the different doseage levels.
#' @param dose_weights Vector specifying weights for the different doses. Default NULL
#' @param prior_list a prior_list object. Default NULL
#' @param sd_posterior tbd. Default NULL
#' @param se_new_trial tbd. Default NULL
#' 
#' @return contr Object of class ‘⁠optContr⁠’. A list containing entries contMat and muMat, and CorrMat. Specified in the Dosefinding package.
#' 
#' @export
getContr <- function (
    
  mods,
  dose_levels,
  dose_weights = NULL,
  prior_list   = NULL,
  sd_posterior = NULL,
  se_new_trial = NULL
  
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
#' @param mods An object of class "Mods" as specified in the Dosefinding package.
#' @param dose_levels vector containing the different dosage levels.
#' @param dose_weights Vector specifying weights for the different doses
#' @param alpha_crit_val significance level. Default set to 0.025.
#' 
#' @return crit_pval multiplicity adjusted critical value on the probability scale.
#' 
#' @export
getCritProb <- function (
    
  mods,
  dose_levels,
  dose_weights,
  alpha_crit_val = 0.025
  
) {
  
  contr <- DoseFinding::optContr(
    models = mods,
    doses  = dose_levels,
    w      = dose_weights)
  
  crit_prob <- stats::pnorm(DoseFinding::critVal(
    corMat      = contr$corMat,
    alpha       = alpha_crit_val,
    df          = 0,
    alternative = "one.sided"))
  
  return (crit_prob)
  
}

#' @title performBayesianMCPMod
#' 
#' @description performs bayesian MCP Test step and modelling.
#' 
#' @param posterior_list a getPosterior object
#' @param contr a getContrMat object, contrast matrix to be used for the testing step.
#' @param crit_prob_adj a getCritProb object
#' @param simple boolean variable, defining whether simplified fit will be applied. Passed to the getModelFits function. Default FALSE.
#' 
#' @return bmcpmod test result as well as modelling result.
#' 
#' @export
performBayesianMCPMod <- function (
    
  posterior_list,
  contr,
  crit_prob_adj,
  simple = FALSE
  
) {
  
  if (inherits(posterior_list,  "postList")) {
    
    posterior_list <- list(posterior_list)
    
  }
  
  if (inherits(contr, "optContr")) {
    
    model_shapes <- colnames(contr$contMat)
    dose_levels  <- as.numeric(rownames(contr$contMat))
    
  } else if (length(contr) == length(posterior_list)) {
    
    model_shapes <- colnames(contr[[1]]$contMat)
    dose_levels  <- as.numeric(rownames(contr[[1]]$contMat))
    
  } else {
    
    stop ("Argument 'contr' must be of type 'optContr'")
    
  }
  
  b_mcp <- performBayesianMCP(
    posterior_list = posterior_list,
    contr          = contr,
    crit_prob_adj  = crit_prob_adj)
  
  fits_list <- lapply(seq_along(posterior_list), function (i) {
    
    if (b_mcp[i, 1]) {
      
      sign_models <- b_mcp[i, -c(1, 2)] > attr(b_mcp, "crit_prob_adj")
      
      model_fits  <- getModelFits(
        models      = model_shapes,
        dose_levels = dose_levels,
        posterior   = posterior_list[[i]],
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
#' @description performs bayesian MCP Test step.
#' 
#' @param posterior_list a getPosterior object
#' @param contr a getContrMat object, contrast matrix to be used for the testing step.
#' @param crit_prob_adj a getCritProb object, specifying the critical value to be used for the testing (on the probability scale)
#' 
#' @return b_mcp test result 
#' 
#' @export
performBayesianMCP <- function(
    
  posterior_list,
  contr,
  crit_prob_adj
  
) {
  
  if (inherits(posterior_list,  "postList")) {
    
    posterior_list <- list(posterior_list)
    
  }
  
  if (inherits(contr, "optContr")) {
    
    b_mcp <- t(sapply(posterior_list, BayesMCPi, contr, crit_prob_adj))
    
  } else {
    
    b_mcp <- t(mapply(BayesMCPi, posterior_list, contr, crit_prob_adj))
    
  }
  
  class(b_mcp)                 <- "BayesianMCP"
  attr(b_mcp, "crit_prob_adj") <- crit_prob_adj
  attr(b_mcp, "ess_avg")       <- ifelse(
    test = is.na(attr(posterior_list[[1]], "ess")),
    yes  = numeric(0),
    no   = rowMeans(sapply(posterior_list, function (posteriors) {
      
      attr(posteriors, "ess")
      
    })))
  
  
  return (b_mcp)
  
}

BayesMCPi <- function (
    
  posterior_i,
  contr,
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
  post_probs   <- apply(contr$contMat, 2, getPostProb, post_combs_i)
  
  res <- c(sign          = ifelse(max(post_probs) > crit_prob_adj, 1, 0),
           crit_prob_adj = crit_prob_adj,
           max_post_prob = max(post_probs),
           post_probs    = post_probs)
  
  return (res)
  
}
