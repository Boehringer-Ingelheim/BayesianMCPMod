addSignificance <- function (

  model_fits,
  sign_models

) {

  sign_models_fits <- sapply(names(model_fits), function (model_name) {
    
    if (model_name != "avgFit") {
      
      indx <- grepl(model_name, names(sign_models))
      
      sum(sign_models[indx]) > 0
      
    } else {
      
      NA
      
    }
    
  })


  model_fits_out <- lapply(seq_along(model_fits), function (i) {

    c(model_fits[[i]], significant = sign_models_fits[i])

  })

  attributes(model_fits_out) <- attributes(model_fits)

  return (model_fits_out)

}

BayesMCPi <- function (

  posterior_i,
  contr,
  crit_prob_adj

) {
  
  post_comps_i <- attr(posterior_i, "posteriorInfo")
  
  ## Pinheiro, José, et al.
  ## "Model‐based dose finding under model uncertainty using general parametric models." 
  ## Statistics in medicine 33.10 (2014): 1646-1661. doi:10.1002/sim.6052
  ## Section 2.2
  
  mu_mat    <- do.call(cbind, post_comps_i$means) # [N candidate models, N posterior components]
  S_list    <- post_comps_i$covMats               # list[[N posterior components]] of cov matrices
                                                  # with [N candidate models, N candidate models]
  
  # Optimal contrast matrix
  contr_mat <- contr$contMat                      # [N dose levels, N candidate models]
  
  ## Contrast estimates for each component
  contr_est_mu  <- t(contr_mat) %*% mu_mat        # [N candidate models, N posterior components]
  contr_est_cov <-
    lapply(S_list, function(S)                    # list[[N posterior components]] of cov matrices
      t(contr_mat) %*% S %*% contr_mat)           # with [N candidate models, N candidate models]
  
  ## Test statistic for each component            # [N candidate models, N posterior components]
  z_m <- contr_est_mu / do.call(cbind, lapply(contr_est_cov, diag))^0.5
  
  ## Fleischer, Frank, et al.
  ## "Bayesian MCPMod."
  ## Pharmaceutical Statistics 21.3 (2022): 654-670. doi:10.1002/pst.2193
  ## Section 2.3
  
  ## Posterior probabilities                      # [N candidate models]
  post_probs <- (stats::pnorm(z_m) %*% unlist(post_comps_i$weights))[, 1]

  res <- c(sign          = ifelse(max(post_probs) > crit_prob_adj, 1, 0),
           crit_prob_adj = crit_prob_adj,
           max_post_prob = max(post_probs),
           post_probs    = post_probs)

  return (res)
  
  # getPostProb <- function (
  #   
  #   contr_j,     # j: dose level
  #   post_combs_i # i: simulation outcome
  #   
  # ) {
  # 
  #   ## Fleischer, Frank, et al. "Bayesian MCPMod."
  #   ## Pharmaceutical Statistics 21.3 (2022): 654-670. doi:10.1002/pst.2193
  #   ## Sections 2.2 and 2.3
  #   
  #   ## Test statistic = sum over all components of
  #   ## posterior weight * normal probability distribution of
  #   ## critical values for doses * estimated mean / sqrt(product of critical values for doses)
  #   
  #   ## Calculation for each component of the posterior
  #   contr_theta   <- apply(post_combs_i$means, 1, `%*%`, contr_j)
  #   contr_var     <- apply(post_combs_i$vars, 1, `%*%`, contr_j^2)
  #   contr_weights <- post_combs_i$weights
  #   
  #   ## P(c_m * theta > 0 | Y = y) for a shape m (and dose j)
  #   post_probs <- sum(contr_weights * stats::pnorm(contr_theta / sqrt(contr_var)))
  #   
  #   return (post_probs)
  #   
  # }
  # 
  # post_combs_i <- getPostCompsI(posterior_i)
  # post_probs   <- apply(contr$contMat, 2, getPostProb, post_combs_i)
  # 
  # res <- c(sign          = ifelse(max(post_probs) > crit_prob_adj, 1, 0),
  #          crit_prob_adj = crit_prob_adj,
  #          max_post_prob = max(post_probs),
  #          post_probs    = post_probs)
  # 
  # return (res)
  
}

#' @title getContr
#'
#' @description This function calculates contrast vectors that are optimal for detecting certain alternatives via applying the function optContr() of the DoseFinding package.
#' Hereby, 4 different options can be distinguished that are automatically executed based on the input that is provided
#' 1) Bayesian approach: If dose_weights and a prior_list are provided an optimized contrasts for the posterior sample size is calculated.
#'    In detail,  in a first step the dose_weights (typically the number of patients per dose group) and the prior information is combined by calculating for
#'    each dose group a posterior effective sample. Based on this posterior effective sample sizes the allocation ratio is derived, which allows for a calculation on
#'    pseudo-optimal contrasts via regular MCPMod are calculated from the
#'    regular MCPMod for these specific weights
#' 2) Frequentist approach: If only dose_weights are provided optimal contrast vectors are calculated from the
#'    regular MCPMod for these specific weights
#' 3) Bayesian approach + re-estimation: If only a cov_posterior (i.e. variability of the posterior distribution) is provided, pseudo-optimal contrasts based on these posterior weights will be calculated
#' 4) Frequentist approach+re-estimation: If only a cov_new_trial (i.e. the estimated variability of a new trial) is provided, optimal contrast vectors are calculated from the
#'    regular MCPMod for this specific covariance matrix. 
#'
#' @param mods An object of class 'Mods' as created by the function 'DoseFinding::Mods()'
#' @param dose_levels Vector containing the different dosage levels.
#' @param dose_weights Vector specifying weights for the different doses. Please note that in case this information is provided together with a prior (i.e. Option 1) is planned these two inputs should be provided on the same scale (e.g. patient numbers).  Default NULL
#' @param prior_list A list of objects of class 'normMix' as created with 'RBesT::mixnorm()'. Only required as input for Option 1. Default NULL
#' @param cov_posterior A covariance matrix with information about the variability of the posterior distribution, only required for Option 3. Default NULL
#' @param cov_new_trial A covariance matrix with information about the observed variability, only required for Option 4. Default NULL
#'
#' @examples
#' dose_levels  <- c(0, 0.5, 2, 4, 8)
#' mods <- DoseFinding::Mods(
#'   linear      = NULL,
#'   emax        = c(0.5, 1.2),
#'   exponential = 2,
#'   doses       = dose_levels,
#'   maxEff      = 6)
#' cov_posterior <- diag(c(2.8, 3, 2.5, 3.5, 4)^2)
#'
#' contr_mat <- getContr(
#'   mods         = mods,
#'   dose_levels  = dose_levels,
#'   cov_posterior = cov_posterior)
#'
#' @return An object of class 'optContr' as provided by the function 'DoseFinding::optContr()'.
#'
#' @export
getContr <- function (
    
    mods,
    dose_levels,
    dose_weights  = NULL,
    prior_list    = NULL,
    cov_posterior = NULL,
    cov_new_trial = NULL
    
) {
  
  checkmate::assert_class(mods, "Mods")
  checkmate::assert_numeric(dose_levels, lower = 0, any.missing = FALSE, len = length(attr(mods, "doses")))

  if (!is.null(dose_weights)) {
    checkmate::assert_numeric(dose_weights, any.missing = FALSE, len = length(attr(mods, "doses")))
  }
  if (!is.null(prior_list)) {
    checkmate::assert_list(prior_list, names = "named", len = length(attr(mods, "doses")), any.missing = FALSE)
  }

  # Determine scenario
  if (!is.null(cov_new_trial) &&
      is.null(dose_weights) && is.null(prior_list) && is.null(cov_posterior)) {

    w <- NULL
    S <- cov_new_trial

  } else if (!is.null(dose_weights) &&
             is.null(cov_new_trial) && is.null(prior_list) && is.null(cov_posterior)) {

    w <- dose_weights
    S <- NULL

  } else if (!is.null(cov_posterior) &&
             is.null(cov_new_trial) && is.null(prior_list) && is.null(dose_weights)) {

    w <- NULL
    S <- cov_posterior

  } else if (!is.null(dose_weights) && !is.null(prior_list) &&
             is.null(cov_new_trial) && is.null(cov_posterior)) {

    w <- dose_weights +
      suppressMessages(round(unlist(lapply(prior_list, RBesT::ess))))
    S <- NULL

  } else {
    stop("Invalid combination of inputs. Valid combinations are:
    1. cov_new_trial only (frequentist with re-estimation)
    2. dose_weights only (frequentist without re-estimation)
    3. cov_posterior only (Bayesian with re-estimation)
    4. dose_weights + prior_list (Bayesian without re-estimation)")
  }

  # Compute contrasts
  if (is.null(w)) {
    contr <- DoseFinding::optContr(
      models = mods,
      doses  = dose_levels,
      S      = S
    )
  } else {
    contr <- DoseFinding::optContr(
      models = mods,
      doses  = dose_levels,
      w      = w
    )
  }
  
  attr(contr, "direction") <- attr(mods, "direction")

  return(contr)
  
}

#' @title getCritProb
#'
#' @description This function calculates multiplicity adjusted critical values. The critical values are calculated in such a way that
#'  when using non-informative priors the actual error level for falsely declaring a significant trial in the Bayesian MCPMod is controlled (by the specified alpha level).
#'  Hereby optimal contrasts of the frequentist MCPMod are applied and two options can be distinguished
#'  1) Frequentist approach: If only dose_weights are provided optimal contrast vectors are calculated from the
#'     regular MCPMod for these specific weights and the corresponding critical value for this set of contrasts is calculated via the critVal() function of the DoseFinding package.
#'  2) Frequentist approach + re-estimation: If only a cov_new_trial (i.e. the covariance matrix of a new trial) is provided, optimal contrast vectors are calculated from the
#'     regular MCPMod for this specific matrix. Here as well the critical value for this set of contrasts is calculated via the critVal() function of the DoseFinding package.
#'
#' @param mods An object of class "Mods" as specified in the DoseFinding package.
#' @param dose_levels Vector containing the different dosage levels.
#' @param dose_weights Vector specifying weights for the different doses, only required for Option i). Default NULL
#' @param cov_new_trial A covariance matrix, only required for Option ii). Default NULL
#' @param alpha_crit_val Significance level. Default set to 0.025.
#'
#' @examples
#' mods <- DoseFinding::Mods(linear      = NULL,
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
  cov_new_trial  = NULL,
  alpha_crit_val = 0.025

) {

  checkmate::assert_class(mods, classes = "Mods")
  checkmate::assert_double(dose_levels, lower = 0, any.missing = FALSE)
  checkmate::assert_double(dose_weights, any.missing = FALSE, len = length(dose_levels), null.ok = TRUE)
  checkmate::assert_double(alpha_crit_val, lower = 0, upper = 1)

  # Get contrast using updated getContr
  contr <- getContr(
    mods          = mods,
    dose_levels   = dose_levels,
    dose_weights  = dose_weights,
    cov_new_trial = cov_new_trial
  )

  # Calculate critical probability
  crit_prob <- stats::pnorm(DoseFinding::critVal(
    corMat      = contr$corMat,
    alpha       = alpha_crit_val,
    df          = 0,
    alternative = "one.sided"
  ))

  return(crit_prob)
  
}

getModelSuccesses <- function (b_mcp) {

  stopifnot(inherits(b_mcp, "BayesianMCP"))

  model_indices <- grepl("post_probs.", colnames(b_mcp))
  model_names   <- colnames(b_mcp)[model_indices] |>
    sub(pattern = "post_probs.", replacement = "", x = _)

  model_successes <- colMeans(b_mcp[, model_indices] > b_mcp[, "crit_prob_adj"])

  names(model_successes) <- model_names

  return (model_successes)

}

#' @title performBayesianMCPMod
#'
#' @description Performs Bayesian MCP Test step and modeling in a combined fashion. See performBayesianMCP() function for MCP Test step and getModelFits() for the modeling step
#'
#' @param posterior_list An object of class 'postList' or a list of 'postList' objects as created by getPosterior() containing information about the (mixture) posterior distribution per dose group
#' @param contr An object of class 'optContr' as created by the getContr() function. It contains the contrast matrix to be used for the testing step.
#' @param crit_prob_adj A getCritProb object, specifying the critical value to be used for the testing (on the probability scale).
#' @param simple Boolean variable, defining whether simplified fit will be applied. Passed to the getModelFits() function. Default FALSE.
#' @param avg_fit Boolean variable, defining whether an average fit (based on generalized AIC weights) should be performed in addition to the individual models. Default TRUE.
#' @param delta A numeric value for the threshold Delta for the MED assessment. If NULL, no MED assessment is performed. Default NULL.
#' @param evidence_level A numeric value between 0 and 1 for the evidence level gamma for the MED assessment. Only required for Bayesian MED assessment, see ?getMED for details. Default NULL.
#' @param med_selection A string, either "avgFit" or "bestFit" based on the lowest gAIC, for the method of MED selection. Default "avgFit".
#' @param n_samples A numerical for the number of bootstrapped samples in case the Bayesian MED assessment is performed. Default 1e3.
#' @param probability_scale A boolean to specify if the trial has a continuous or a binary outcome. Setting to TRUE will transform calculations from the logit scale to the probability scale, which can be desirable for a binary outcome. Default FALSE.
#' @examples
#' mods <- DoseFinding::Mods(linear      = NULL,
#'                           emax        = c(0.5, 1.2),
#'                           exponential = 2,
#'                           doses       = c(0, 0.5, 2,4, 8))
#' dose_levels  <- c(0, 0.5, 2, 4, 8)
#' sd_posterior <- c(2.8, 3, 2.5, 3.5, 4)
#' contr_mat <- getContr(
#'   mods          = mods,
#'   dose_levels   = dose_levels,
#'   cov_posterior = diag(sd_posterior)^2)
#' critVal <- getCritProb(
#'   mods           = mods,
#'   dose_weights   = c(50, 50, 50, 50, 50), #reflecting the planned sample size
#'   dose_levels    = dose_levels,
#'   alpha_crit_val = 0.6) # unreasonable alpha chosen for this example, rather choose 0.05
#' prior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 5), sigma = 2),
#'                    DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 1, s = 12), sigma = 2),
#'                    DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.2, s = 11), sigma = 2) ,
#'                    DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.3, s = 11), sigma = 2) ,
#'                    DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 2, s = 13), sigma = 2))
#' mu <- c(0, 1, 1.5, 2, 2.5)
#' S_hat <- diag(c(5, 4, 6, 7, 8)^2)
#' posterior_list <- getPosterior(
#'   prior_list = prior_list,
#'   mu_hat     = mu,
#'   S_hat      = S_hat,
#'   calc_ess   = TRUE)
#'
#' performBayesianMCPMod(posterior_list = posterior_list,
#'                       contr          = contr_mat,
#'                       crit_prob_adj  = critVal,
#'                       simple         = FALSE,
#'                       delta          = 1.1)
#'
#' @return Bayesian MCP test result as well as modeling result.
#'
#' @export
performBayesianMCPMod <- function (

  posterior_list,
  contr,
  crit_prob_adj,
  simple            = FALSE,
  avg_fit           = TRUE,

  delta             = NULL,
  evidence_level    = NULL,
  med_selection     = c("avgFit", "bestFit"),
  n_samples         = 1e3,
  probability_scale = FALSE

) {

  if (inherits(posterior_list,  "postList")) {
    posterior_list <- list(posterior_list)
  }
  checkmate::assert_list(posterior_list, types = "postList")
  
  checkmate::assert_class(contr, "optContr")
  checkmate::assert_class(crit_prob_adj, "numeric")
  checkmate::assert_logical(simple)
  checkmate::assert_logical(avg_fit)

  checkmate::assert_double(delta, null.ok = TRUE)
  checkmate::assert_double(evidence_level, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_character(med_selection)
  checkmate::assert_double(n_samples, lower = 1)

  med_selection <- match.arg(med_selection, choices = c("avgFit", "bestFit"))
  get_med       <- ifelse(is.null(delta), FALSE, TRUE)
  
  checkmate::assert_flag(probability_scale)

  if (!is.null(evidence_level)) {
    stopifnot("delta must not be NULL if evidence_level is not NULL" = !is.null(delta))
  }

  if (inherits(contr, "optContr")) {

    model_names <- colnames(contr$contMat)
    dose_levels <- as.numeric(rownames(contr$contMat))

  } else if (length(contr) == length(posterior_list)) {

    model_names <- colnames(contr[[1]]$contMat)
    dose_levels <- as.numeric(rownames(contr[[1]]$contMat))

  } else {

    stop ("Argument 'contr' must be of type 'optContr'")

  }
  
  if (is.null(attr(contr, "direction"))) {
    
    attr(model_names, "direction") <- "increasing"
    message("'contr' does not have attribute 'direction' and is assumed to be 'increasing'")
    
  } else {
    
    attr(model_names, "direction") <- attr(contr, "direction")
    
  }

  b_mcp <- performBayesianMCP(
    posterior_list = posterior_list,
    contr          = contr,
    crit_prob_adj  = crit_prob_adj)

  b_mod <- performBayesianMod(
    b_mcp             = b_mcp,
    posterior_list    = posterior_list,
    model_names       = model_names,
    dose_levels       = dose_levels,
    simple            = simple,
    avg_fit           = avg_fit,
    probability_scale = probability_scale)

  b_mcp_mod <- list(BayesianMCP = b_mcp, Mod = b_mod)

  if (get_med) {

    med_info <- t(sapply(b_mcp_mod$Mod, function (model_fits) {

      if (!is.null(model_fits)) { # model shape is significant

        if (!is.null(evidence_level)) { # && abs(evidence_level - 0.5) < 1e-9 # in case MED assessment with 0.5 should result in using the modelFit
          
          bs_quantile <- ifelse(test = attr(model_fits, "direction") == "increasing",
                                yes  = 1 - evidence_level,
                                no   = evidence_level)

          bs_quantiles <- getBootstrapQuantiles(
            model_fits        = model_fits,
            quantiles         = bs_quantile,
            n_samples         = n_samples,
            probability_scale = probability_scale)

          med_info <- getMED(
            delta             = delta,
            evidence_level    = evidence_level,
            bs_quantiles      = bs_quantiles,
            probability_scale = probability_scale)

        } else {

          med_info <- getMED(
            delta             = delta,
            model_fits        = model_fits,
            probability_scale = probability_scale)

        }

        if (med_selection == "avgFit") {

          return (med_info[, "avgFit"])

        } else if (med_selection == "bestFit") {

          model_min <- which.min(sapply(model_fits, function (x) x$gAIC))

          return (med_info[, model_min])

        }

      } else {

        return (c(med_reached  = 0, med = NA))

      }

    }))

    attr(b_mcp_mod, "MED")          <- med_info
    attr(b_mcp_mod, "MEDSelection") <- med_selection

  }

  class(b_mcp_mod) <- "BayesianMCPMod"

  return (b_mcp_mod)

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
#'                           emax        = c(0.5, 1.2),
#'                           exponential = 2,
#'                           doses       = c(0, 0.5, 2,4, 8))
#' dose_levels  <- c(0, 0.5, 2, 4, 8)
#' sd_posterior <- c(2.8,3,2.5,3.5,4)
#' contr_mat <- getContr(
#'   mods          = mods,
#'   dose_levels   = dose_levels,
#'   cov_posterior = diag(sd_posterior)^2)
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
#' S_hat <- diag(c(5, 4, 6, 7, 8)^2)
#' posterior_list <- getPosterior(
#'   prior_list = prior_list,
#'   mu_hat     = mu,
#'   S_hat      = S_hat,
#'   calc_ess   = TRUE)
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
  
  if (inherits(posterior_list,  "postList")) {
    
    posterior_list <- list(posterior_list)
    
  }
  checkmate::assert_list(posterior_list, types = "postList")
  
  checkmate::assert(
    checkmate::check_class(contr, "optContr"),
    checkmate::check_list(contr, types = "optContr"),
    combine = "or")
  
  checkmate::assert_class(crit_prob_adj, "numeric")
  checkmate::assert_numeric(crit_prob_adj, lower = 0, upper = Inf)

  ## What case is covered by the IF clause?
  ## -> if contrasts are re-evaluated, contr becomes a list of contrasts
  if (inherits(contr, "optContr")) {

    b_mcp <- t(sapply(posterior_list, BayesMCPi, contr, crit_prob_adj))

  } else {

    b_mcp <- t(mapply(BayesMCPi, posterior_list, contr, crit_prob_adj))

  }

  class(b_mcp)               <- "BayesianMCP"
  attr(b_mcp, "critProbAdj") <- crit_prob_adj
  attr(b_mcp, "successRate") <- mean(b_mcp[, 1])
  attr(b_mcp, "essAvg")      <- ifelse(
    test = is.na(attr(posterior_list[[1]], "ess")),
    yes  = numeric(0),
    no   = rowMeans(sapply(posterior_list,
                           function (posteriors) attr(posteriors, "ess"))))

  return (b_mcp)

}

performBayesianMod <- function (

  b_mcp,
  posterior_list,
  model_names,
  dose_levels,
  simple,
  avg_fit,
  probability_scale = FALSE

) {

  fits_list <- optPar_lapply(seq_along(posterior_list), function (i) {

    if (b_mcp[i, 1]) {

      model_fits <- getModelFits(
        models            = model_names,
        dose_levels       = dose_levels,
        posterior         = posterior_list[[i]],
        avg_fit           = avg_fit,
        simple            = simple,
        probability_scale = probability_scale)

      sign_models <- b_mcp[i, -c(1, 2)] > attr(b_mcp, "critProbAdj")
      model_fits  <- addSignificance(model_fits, sign_models)

    } else {

      NULL

    }

  })
  
  # requireNamespace("DoseFinding", quietly = TRUE)
  # requireNamespace("RBesT",       quietly = TRUE)
  # requireNamespace("nloptr",      quietly = TRUE)
  # 
  # 
  # fits_list <- parallel::mclapply(seq_along(posterior_list), function (i) {
  #   
  #   if (b_mcp[i, 1]) {
  #     
  #     model_fits <- getModelFits(
  #       models            = model_names,
  #       dose_levels       = dose_levels,
  #       posterior         = posterior_list[[i]],
  #       avg_fit           = avg_fit,
  #       simple            = simple,
  #       probability_scale = probability_scale)
  #     
  #     sign_models <- b_mcp[i, -c(1, 2)] > attr(b_mcp, "critProbAdj")
  #     model_fits  <- addSignificance(model_fits, sign_models)
  #     
  #   } else {
  #     
  #     NULL
  #     
  #   }
  #   
  # },
  # mc.cores       = 32L,
  # mc.preschedule = TRUE,
  # mc.set.seed    = FALSE)

}
