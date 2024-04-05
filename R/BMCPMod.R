#' @title assessDesign
#'.
#' @description This function performs simulation based trial design evaluations for a set of specified dose-response models
#'
#' @param n_patients Vector specifying the planned number of patients per dose group
#' @param mods An object of class "Mods" as specified in the DoseFinding package.
#' @param prior_list A prior_list object specifying the utilized prior for the different dose groups 
#' @param sd A positive value, specification of assumed sd 
#' @param n_sim Number of simulations to be performed
#' @param alpha_crit_val (Un-adjusted) Critical value to be used for the MCP testing step. Passed to the getCritProb() function for the calculation of adjusted critical values (on the probability scale). Default is 0.05.
#' @param simple Boolean variable defining whether simplified fit will be applied. Passed to the getModelFits function. Default FALSE.
#' @param reestimate Boolean variable defining whether critical value should be calculated with re-estimated contrasts (see getCritProb function for more details). Default FALSE
#' @param contr An object of class 'optContr' as created by the getContr() function. Allows specification of a fixed contrasts matrix. Default NULL
#' @param dr_means A vector, allows specification of  individual (not model based) assumed effects per dose group. Default NULL
#' 
#' @return Returns success probabilities for the different assumed dose-response shapes, attributes also includes information around average success rate (across all assumed models) and prior Effective sample size
#' 
#' @examples
#' if (interactive()) { # takes typically > 5 seconds
#' 
#' mods <- DoseFinding::Mods(linear      = NULL,
#'                           linlog      = NULL,
#'                           emax        = c(0.5, 1.2),
#'                           exponential = 2,
#'                           doses       = c(0, 0.5, 2,4, 8),
#'                           maxEff      = 6)
#' sd <- 12
#' prior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 12), sigma = 2),
#'                    DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 1, s = 12), sigma = 2),
#'                    DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.2, s = 11), sigma = 2) ,  
#'                    DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.3, s = 11), sigma = 2) ,
#'                    DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 2, s = 13), sigma = 2))
#' n_patients <- c(40, 60, 60, 60, 60)
#' 
#' success_probabilities <- assessDesign(
#'   n_patients  = n_patients,
#'   mods        = mods,
#'   prior_list  = prior_list,
#'   sd          = sd,
#'   n_sim       = 1e2) # speed up exammple run time
#' 
#' success_probabilities
#' 
#' }
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
  
  checkmate::assert_vector(n_patients, len = length(attr(mods, "doses")), any.missing = FALSE)
  checkmate::check_class(mods, classes = "Mods")
  checkmate::check_list(prior_list, names = "named", len = length(attr(mods, "doses")), any.missing = FALSE)
  # sensitive to how DoseFinding labels their attributes for "Mods" class
  checkmate::check_double(n_sim, lower = 1, upper = Inf)
  checkmate::check_double(alpha_crit_val, lower = 0, upper = 1)
  checkmate::check_logical(simple)
  # TODO: check that prior_list has 'sd_tot' attribute, and that it's numeric # this is not applicable at the moment
  
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
  
  avg_success_rate <- mean(sapply(eval_design, function (bmcpmod) {
    attr(bmcpmod$BayesianMCP, "successRate")
  }))
  
  names(eval_design) <- model_names
  
  attr(eval_design, "avgSuccessRate") <- avg_success_rate
  attr(eval_design, "placEff")        <- ifelse(test = is.null(dr_means),
                                                yes  = attr(mods, "placEff"),
                                                no   = dr_means[1])
  attr(eval_design, "maxEff")         <- ifelse(test = is.null(dr_means),
                                                yes  = attr(mods, "maxEff"),
                                                no   = diff(range(dr_means)))
  attr(eval_design, "sampleSize")     <- n_patients
  attr(eval_design, "priorESS")       <- round(getESS(prior_list), 1)
  
  return (eval_design)
  
}

#' @title getContr
#' 
#' @description This function calculates contrast vectors that are optimal for detecting certain alternatives via applying the function optContr() of the DoseFinding package.
#' Hereby 4 different options can be distinguished that are automatically executed based on the input that is provided
#' 1) Bayesian approach: If dose_weights and a prior_list are provided an optimized contrasts for the posterior sample size is calculated. 
#'    In detail,  in a first step the dose_weights (typically the number of patients per dose group) and the prior information is combined by calculating for
#'    each dose group a posterior effective sample. Based on this posterior effective sample sizes the allocation ratio is derived, which allows for a calculation on
#'    pseudo-optimal contrasts via regular MCPMod are calculated from the
#'    regular MCPMod for these specific weights 
#' 2) Frequentist approach: If only dose_weights are provided optimal contrast vectors are calculated from the
#'    regular MCPMod for these specific weights
#' 3) Bayesian approach + re-estimation: If only a sd_posterior (i.e. variability of the posterior distribution) is provided, pseudo-optimal contrasts based on these posterior weights will be calculated
#' 4) Frequentist approach+re-estimation:If only a se_new_trial (i.e. the estimated variability per dose group of a new trial) is provided, optimal contrast vectors are calculated from the
#'    regular MCPMod for this specific vector of standard errors. For the actual evaluation this vector of standard errors is translated into a (diagonal) matrix of variances 
#' 
#' @param mods An object of class 'Mods' as created by the function 'DoseFinding::Mods()'
#' @param dose_levels Vector containing the different dosage levels.
#' @param dose_weights Vector specifying weights for the different doses. Please note that in case this information is provided together with a prior (i.e. Option 1) is planned these two inputs should be provided on the same scale (e.g. patient numbers).  Default NULL
#' @param prior_list A list of objects of class 'normMix' as created with 'RBesT::mixnorm()'. Only required as input for Option 1. Default NULL
#' @param sd_posterior A vector of positive values with information about the variability of the posterior distribution, only required for Option 3. Default NULL
#' @param se_new_trial A vector of positive values with information about the observed variability, only required for Option 4. Default NULL

#' 
#' @examples
#' dose_levels  <- c(0, 0.5, 2, 4, 8)
#' mods <- DoseFinding::Mods(
#'   linear      = NULL,
#'   linlog      = NULL,
#'   emax        = c(0.5, 1.2),
#'   exponential = 2,
#'   doses       = dose_levels,
#'   maxEff      = 6)
#' sd_posterior <- c(2.8, 3, 2.5, 3.5, 4)
#' 
#' contr_mat <- getContr(
#'   mods         = mods,
#'   dose_levels  = dose_levels,
#'   sd_posterior = sd_posterior) 
#' 
#' @return An object of class 'optContr' as provided by the function 'DoseFinding::optContr()'.
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
  
  checkmate::check_class(mods, classes = "Mods")
  checkmate::check_double(dose_levels, lower = 0, any.missing = FALSE, len = length(attr(mods, "doses")))
  checkmate::check_double(dose_weights, any.missing = FALSE, len = length(attr(mods, "doses")))
  checkmate::check_list(prior_list, names = "named", len = length(attr(mods, "doses")), any.missing = FALSE)
  
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
#' @description This function calculates multiplicity adjusted critical values. The critical values are calculated in such a way that
#'  when using non-informative priors the actual error level for falsely declaring a significant trial in the Bayesian MCPMod is controlled (by the specified alpha level). 
#'  Hereby optimal contrasts of the frequentist MCPMod are applied and two options can be distinguished
#'  1) Frequentist approach: If only dose_weights are provided optimal contrast vectors are calculated from the
#'     regular MCPMod for these specific weights and the corresponding critical value for this set of contrasts is calculated via the critVal() function of the DoseFinding package.
#'  2) Frequentist approach + re-estimation: If only a se_new_trial (i.e. the estimated variability per dose group of a new trial) is provided, optimal contrast vectors are calculated from the
#'     regular MCPMod for this specific vector of standard errors. Here as well the critical value for this set of contrasts is calculated via the critVal() function of the DoseFinding package.
#' 
#' @param mods An object of class "Mods" as specified in the DoseFinding package.
#' @param dose_levels Vector containing the different dosage levels.
#' @param dose_weights Vector specifying weights for the different doses, only required for Option i). Default NULL
#' @param se_new_trial A vector of positive values, only required for Option ii). Default NULL
#' @param alpha_crit_val Significance level. Default set to 0.025.
#'
#' @examples
#' mods <- DoseFinding::Mods(linear      = NULL,
#'                           linlog      = NULL,
#'                           emax        = c(0.5, 1.2),
#'                           exponential = 2,
#'                           doses       = c(0, 0.5, 2,4, 8))
#' dose_levels <- c(0, 0.5, 2, 4, 8)
#' critVal <- getCritProb(
#'   mods           = mods,
#'   dose_weights   = c(50,50,50,50,50), #reflecting the planned sample size
#'   dose_levels    = dose_levels,
#'   alpha_crit_val = 0.05) 
#' @return Multiplicity adjusted critical value on the probability scale.
#' 
#' @export
getCritProb <- function (
    
  mods,
  dose_levels,
  dose_weights   = NULL,
  se_new_trial   = NULL,
  alpha_crit_val = 0.025
  
) {
  
  checkmate::check_class(mods, classes = "Mods")
  checkmate::check_double(dose_levels, lower = 0, any.missing = FALSE, len = length(dose_weights))
  checkmate::check_double(dose_weights, any.missing = FALSE, len = length(dose_levels))
  checkmate::check_double(alpha_crit_val, lower = 0, upper = 1)
  
  contr <- getContr(mods           = mods,
                    dose_levels    = dose_levels ,
                    dose_weights   = dose_weights,
                    se_new_trial   = se_new_trial)
  
  crit_prob <- stats::pnorm(DoseFinding::critVal(
    corMat      = contr$corMat,
    alpha       = alpha_crit_val,
    df          = 0,
    alternative = "one.sided"))
  
  return (crit_prob)
  
}

#' @title performBayesianMCPMod
#' 
#' @description Performs Bayesian MCP Test step and modeling in a combined fashion. See performBayesianMCP() function for MCP Test step and getModelFits() for the modeling step
#' 
#' @param posterior_list An object of class 'postList' as created by getPosterior() containing information about the (mixture) posterior distribution per dose group
#' @param contr An object of class 'optContr' as created by the getContr() function. It contains the contrast matrix to be used for the testing step.
#' @param crit_prob_adj A getCritProb object, specifying the critical value to be used for the testing (on the probability scale).
#' @param simple Boolean variable, defining whether simplified fit will be applied. Passed to the getModelFits() function. Default FALSE.
#' @examples
#' mods <- DoseFinding::Mods(linear      = NULL,
#'                           linlog      = NULL,
#'                           emax        = c(0.5, 1.2),
#'                           exponential = 2,
#'                           doses       = c(0, 0.5, 2,4, 8))
#' dose_levels  <- c(0, 0.5, 2, 4, 8)
#' sd_posterior <- c(2.8, 3, 2.5, 3.5, 4)
#' contr_mat <- getContr(
#'   mods         = mods,
#'   dose_levels  = dose_levels,
#'   sd_posterior = sd_posterior)
#' critVal <- getCritProb(
#'   mods           = mods,
#'   dose_weights   = c(50, 50, 50, 50, 50), #reflecting the planned sample size
#'   dose_levels    = dose_levels,
#'   alpha_crit_val = 0.05)
#' prior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 5), sigma = 2),
#'                    DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 1, s = 12), sigma = 2),
#'                    DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.2, s = 11), sigma = 2) ,  
#'                    DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.3, s = 11), sigma = 2) ,
#'                    DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 2, s = 13), sigma = 2))
#' mu <- c(0, 1, 1.5, 2, 2.5)
#' se <- c(5, 4, 6, 7, 8)
#' posterior_list <- getPosterior(
#'   prior_list = prior_list,
#'   mu_hat     = mu,
#'   se_hat     = se)
#' performBayesianMCPMod(posterior_list = posterior_list,
#'                       contr          = contr_mat,
#'                       crit_prob_adj  = critVal,
#'                       simple         = FALSE)
#' 
#' @return Bayesian MCP test result as well as modelling result.
#' 
#' @export
performBayesianMCPMod <- function (
    
  posterior_list,
  contr,
  crit_prob_adj,
  simple = FALSE
  
) {
  
  checkmate::check_class(posterior_list, "postList")
  checkmate::check_class(contr, "optContr")
  checkmate::check_class(crit_prob_adj, "numeric")
  checkmate::check_logical(simple)
  
  if (inherits(posterior_list,  "postList")) {
    
    posterior_list <- list(posterior_list)
    
  }
  
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
#' @description Performs Bayesian MCP Test step, as described in Fleischer et al. (2022).
#' Tests for a dose-response effect using a model-based multiple contrast test based on the (provided) posterior distribution. In particular for every dose-response candidate the posterior probability is calculated that the contrast is bigger than 0 (based on the posterior distribution of the dose groups).
#' In order to obtain significant test decision we consider the maximum of the posterior probabilities across the different models. This maximum is compared with a (multiplicity adjusted) critical value (on the probability scale).
#' @references Fleischer F, Bossert S, Deng Q, Loley C, Gierse J. 2022. Bayesian MCPMod. Pharmaceutical Statistics. 21(3): 654-670. doi:10.1002/pst.2193
#' @param posterior_list An object derived with getPosterior with information about the (mixture) posterior distribution per dose group 
#' @param contr An object of class 'optContr' as created by the getContr() function. It contains the contrast matrix to be used for the testing step.
#' @param crit_prob_adj A getCritProb object, specifying the critical value to be used for the testing (on the probability scale)
#' 
#' @examples
#' mods <- DoseFinding::Mods(linear      = NULL,
#'                           linlog      = NULL,
#'                           emax        = c(0.5, 1.2),
#'                           exponential = 2,
#'                           doses       = c(0, 0.5, 2,4, 8))
#' dose_levels  <- c(0, 0.5, 2, 4, 8)
#' sd_posterior <- c(2.8,3,2.5,3.5,4)
#' contr_mat <- getContr(
#'   mods         = mods,
#'   dose_levels  = dose_levels,
#'   sd_posterior = sd_posterior)
#' critVal <- getCritProb(
#'   mods           = mods,
#'   dose_weights   = c(50, 50, 50, 50, 50), #reflecting the planned sample size
#'   dose_levels    = dose_levels,
#'   alpha_crit_val = 0.05)
#' prior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 5), sigma = 2),
#'                    DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 1, s = 12), sigma = 2),
#'                    DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.2, s = 11), sigma = 2) ,  
#'                    DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.3, s = 11), sigma = 2) ,
#'                    DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 2, s = 13), sigma = 2))
#' mu <- c(0, 1, 1.5, 2, 2.5)
#' se <- c(5, 4, 6, 7, 8)
#' posterior_list <- getPosterior(
#'   prior_list = prior_list,
#'   mu_hat     = mu,
#'   se_hat     = se)
#'   
#' performBayesianMCP(posterior_list = posterior_list,
#'                    contr          = contr_mat,
#'                    crit_prob_adj  = critVal)
#' 
#' @return Bayesian MCP test result, with information about p-values for the individual dose-response shapes and overall significance    
#' 
#' @export
performBayesianMCP <- function(
    
  posterior_list,
  contr,
  crit_prob_adj
  
) {
  
  checkmate::check_class(posterior_list, "postList")
  checkmate::check_class(contr, "optContr")
  checkmate::check_class(crit_prob_adj, "numeric")
  checkmate::check_numeric(crit_prob_adj, lower = 0, upper = Inf)
  
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
  attr(b_mcp, "successRate")   <- mean(b_mcp[, 1])
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