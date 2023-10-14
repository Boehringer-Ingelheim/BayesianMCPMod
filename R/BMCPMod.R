#' @title assessDesign
#' 
#' @param n_patients tbd
#' @param dose_levels tbd
#' @param sd tbd
#' @param mods tbd
#' @param prior_list tbd
#' @param n_sim tbd
#' @param alpha_crit_val tbd
#' @param simple tbd
#' 
#' @export
assessDesign <- function (
    
  n_patients,
  dose_levels,
  sd,
  mods,
  prior_list,
  
  n_sim          = 1e3,
  alpha_crit_val = 0.05,
  simple         = TRUE
  
) {
  
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
    
    contr_mat <- DoseFinding::optContr(
      models = mods,
      doses  = dose_levels,
      w      = n_patients)
    
    crit_pval <- pnorm(DoseFinding:::critVal(
      corMat      = contr_mat$corMat,
      alpha       = alpha_crit_val,
      df          = 0,
      alternative = "one.sided"))
    
    ess_prior <- suppressMessages(round(unlist(lapply(prior_list, RBesT::ess))))
    contr_mat_prior <- DoseFinding::optContr(
      models = mods,
      doses  = dose_levels,
      w      = n_patients + ess_prior)
    
    b_mcp_mod <- performBayesianMCPMod(
      posteriors_list = posterior_list,
      contr_mat       = contr_mat_prior,
      crit_prob       = crit_pval,
      simple          = simple)
    
  })
  
  names(eval_design) <- model_names
  
  return (eval_design)
  
}

#' @title performBayesianMCPMod
#' 
#' @param posteriors_list tbd
#' @param contr_mat tbd
#' @param crit_prob tbd
#' @param simple tbd
#' 
#' @export
performBayesianMCPMod <- function (
    
  posteriors_list,
  contr_mat,
  crit_prob,
  simple = FALSE
  
) {
  
  if (class(posteriors_list) == "postList") {
    
    posteriors_list <- list(posteriors_list)
    
  }
  
  b_mcp <- performBayesianMCP(
    posteriors_list = posteriors_list,
    contr_mat       = contr_mat,
    crit_prob       = crit_prob)
  
  model_shapes <- colnames(contr_mat$contMat)
  dose_levels  <- as.numeric(rownames(contr_mat$contMat))
  
  fits_list <- lapply(seq_along(posteriors_list), function (i) {
    
    if (b_mcp[i, 1]) {
      
      sign_models <- b_mcp[i, -c(1, 2)] > attr(b_mcp, "crit_prob")
      
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

#' @title BayesianMCP
#' 
#' @param posteriors_list tbd
#' @param contr_mat tbd
#' @param crit_prob tbd
#' 
#' @export
performBayesianMCP <- function(
    
  posteriors_list,
  contr_mat,
  crit_prob
  
) {
  
  if (class(posteriors_list) == "postList") {
    
    posteriors_list <- list(posteriors_list)
    
  }
  
  b_mcp <- t(sapply(posteriors_list, BayesMCPi, contr_mat, crit_prob))
  
  attr(b_mcp, "crit_prob") <- crit_prob
  class(b_mcp)             <- "BayesianMCP"
  
  return (b_mcp)
  
}

BayesMCPi <- function (
    
  posterior_i,
  contr_mat,
  crit_prob
  
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
  
  res <- c(sign       = ifelse(max(post_probs) > crit_prob, 1, 0),
           p_val      = max(post_probs),
           post_probs = post_probs)
  
  return (res)
  
}
