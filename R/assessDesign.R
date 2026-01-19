#' @title assessDesign
#'
#' @description This function performs simulation based trial design evaluations for a set of specified dose-response models
#'
#' @param n_patients Vector specifying the planned number of patients per dose group. A minimum of 2 patients are required in each group.
#' @param mods An object of class "Mods" as specified in the DoseFinding package.
#' @param prior_list A prior_list object specifying the utilized prior for the different dose groups
#' @param sd A positive value, specification of assumed sd. Not required if ´data_sim´ or ´estimates_sim´ is provided. Default NULL
#' @param contr An object of class 'optContr' as created by the getContr() function. Allows specification of a fixed contrasts matrix. Default NULL.
#' @param dr_means A vector, allows specification of individual (not model based) assumed effects per dose group. Default NULL.
#' @param data_sim An optional data frame for custom simulated data. Must follow the data structure as provided by ´simulateData()´. Default NULL.
#' @param estimates_sim An optional named list of 1) list of vectors for the estimated means per dose group (estimates_sim$mu_hats) and 2) of list of matrices for the covariance matrices specifying the (estimated) variabilities (estimates_sim$S_hats). Dimensions of entries must match the number of dose levels. Default NULL.
#' @param n_sim Number of simulations to be performed
#' @param alpha_crit_val (Un-adjusted) Critical value to be used for the MCP testing step. Passed to the getCritProb() function for the calculation of adjusted critical values (on the probability scale). Default 0.05.
#' @param modeling Boolean variable defining whether the Mod part of Bayesian MCP-Mod will be performed in the assessment. More heavy on resources. Default FALSE.
#' @param simple Boolean variable defining whether simplified fit will be applied. Passed to the getModelFits function. Default FALSE.
#' @param avg_fit Boolean variable, defining whether an average fit (based on generalized AIC weights) should be performed in addition to the individual models. Default TRUE.
#' @param reestimate Boolean variable defining whether critical value should be calculated with re-estimated contrasts (see getCritProb function for more details). Default FALSE.
#' @param delta A numeric value for the threshold Delta for the MED assessment. If NULL, no MED assessment is performed. Default NULL.
#' @param evidence_level A numeric value between 0 and 1 for the evidence level gamma for the MED assessment. Only required for Bayesian MED assessment, see ?getMED for details. Default NULL.
#' @param med_selection A string, either "avgFit" or "bestFit", for the method of MED selection. Default "avgFit".
#'
#' @return Returns success probabilities for the different assumed dose-response shapes, attributes also includes information around average success rate (across all assumed models) and prior Effective sample size.
#'
#' @examples
#' mods <- DoseFinding::Mods(linear      = NULL,
#'                           emax        = c(0.5, 1.2),
#'                           exponential = 2,
#'                           betaMod     = c(1, 1),
#'                           doses       = c(0, 0.5, 2,4, 8),
#'                           maxEff      = 6)
#'                           
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
#'   n_sim       = 1e2) # speed up example run time
#'
#' success_probabilities
#' 
#' # Analysis with custom dose response relationship
#' custom_dr_means <- c(1, 2, 3, 4, 5)
#' 
#' success_probs_custom_dr <- assessDesign(
#'   n_patients  = n_patients,
#'   mods        = mods,
#'   prior_list  = prior_list,
#'   dr_means    = custom_dr_means,
#'   sd          = sd,
#'   n_sim       = 1e2) # speed up example run time
#'
#' success_probs_custom_dr
#' 
#' # Analysis with custom estimates for means and variabilies
#' # No simulated data, only simulated model estimates
#' estimates_sim <- list(mu_hats = replicate(100, list(c(1, 2, 3, 4, 5) + rnorm(5, 0, 1))),
#'                       S_hats  = list(diag(1, 5)))
#' 
#' success_probs_custom_est <- assessDesign(
#'   n_patients    = n_patients,
#'   mods          = mods,
#'   prior_list    = prior_list,
#'   estimates_sim = estimates_sim)
#'
#' success_probs_custom_est
#'
#' if (interactive()) { # takes typically > 5 seconds
#'
#' # with MED estimation without bootstrapping
#' # see ?getMED for details
#'
#' success_probabilities <- assessDesign(
#'   n_patients     = n_patients,
#'   mods           = mods,
#'   prior_list     = prior_list,
#'   sd             = sd,
#'   modeling       = TRUE,
#'   n_sim          = 10, # speed up example run time
#'   delta          = 7)
#'
#'   success_probabilities
#'
#' # with MED estimation with bootstrapping
#'
#' success_probabilities <- assessDesign(
#'   n_patients     = n_patients,
#'   mods           = mods,
#'   prior_list     = prior_list,
#'   sd             = sd,
#'   modeling       = TRUE,
#'   n_sim          = 10, # speed up example run time
#'   delta          = 7,
#'   evidence_level = 0.8)
#'
#'   success_probabilities
#'
#' }
#'
#' @export
assessDesign <- function (
    
  n_patients,
  mods,
  prior_list,
  
  sd             = NULL,
  
  contr          = NULL,
  dr_means       = NULL,
  
  data_sim       = NULL,
  estimates_sim  = NULL,
  
  n_sim          = 1e3,
  alpha_crit_val = 0.05,
  modeling       = FALSE,
  simple         = TRUE,
  avg_fit        = TRUE,
  reestimate     = FALSE,
  
  delta          = NULL,
  evidence_level = NULL,
  med_selection  = c("avgFit", "bestFit")
  
) {
  
  checkmate::assert_vector(n_patients, len = length(attr(mods, "doses")), any.missing = FALSE)
  checkmate::assert_double(n_patients, lower = 2, upper = Inf)
  checkmate::assert_class(mods, classes = "Mods")
  checkmate::assert_list(prior_list, names = "named", len = length(attr(mods, "doses")), any.missing = FALSE)
  checkmate::assert_list(estimates_sim, null.ok = TRUE)
  # sensitive to how DoseFinding labels their attributes for "Mods" class
  checkmate::assert_integerish(n_sim, lower = 1, upper = Inf)
  checkmate::assert_double(alpha_crit_val, lower = 0, upper = 1)
  checkmate::assert_logical(modeling)
  
  # TODO: check that prior_list has 'sd_tot' attribute, and that it's numeric # this is not applicable at the moment
  
  stopifnot("Cannot specify both 'data_sim' and 'estimates_sim' simultaneously" = 
              is.null(data_sim)  & is.null(estimates_sim) |
              !is.null(data_sim) & is.null(estimates_sim) |
              is.null(data_sim)  & !is.null(estimates_sim))
  
  modeling <- ifelse(!is.null(delta), TRUE, modeling)
  
  dose_levels <- attr(mods, "doses")
  
  if (!is.null(estimates_sim)) {
    
    stopifnot("'estimates_sim' must be a named list of lenght 2 with list items called 'mu_hats' and 'S_hats'" =
                length(estimates_sim) ==  2L &
                all(names(estimates_sim) %in% c("mu_hats", "S_hats")))
    stopifnot("The length of 'estimates_sim$S_hats' must match 'estimates_sim$mu_hats' or 1." =
                length(estimates_sim$S_hats) == length(estimates_sim$mu_hats) |
                length(estimates_sim$S_hats) == 1L)
    stopifnot("The length of each entry of 'estimates_sim$mu_hats' must match the number of dose levels." =
                length(estimates_sim$mu_hats[[1]]) == length(dose_levels))
    stopifnot("The dimensions of each entry of 'estimates_sim$S_hats' must match the number of dose levels." =
                ncol(estimates_sim$S_hats[[1]]) == length(dose_levels))
    
  } else if (!is.null(data_sim)) {
    
    ## lazily simulate data anyway ...
    # TODO: this approach requires people to simulate data for each provided model in mods
    data <- simulateData(
      n_patients  = n_patients,
      dose_levels = dose_levels,
      sd          = 1,
      mods        = mods,
      n_sim       = nrow(data_sim), # adjust the number of simulations to data_sim
      dr_means    = dr_means)
    
    ## ... and then check if formats of data_sim and data match
    stopifnot("'data_sim' must follow the data structure as provided by simulateData()" =
                data_sim[, 1:3] == data[, 1:3])
    
    data <- data_sim
    
  } else {
    
    if (is.null(sd)) stop ("Must provide 'sd' argument for símulation.")
    
    data <- simulateData(
      n_patients  = n_patients,
      dose_levels = dose_levels,
      sd          = sd,
      mods        = mods,
      n_sim       = n_sim,
      dr_means    = dr_means)
    
  }
  
  crit_prob_adj <- getCritProb(
    mods           = mods,
    dose_levels    = dose_levels,
    dose_weights   = n_patients,
    alpha_crit_val = alpha_crit_val)
  
  if (!reestimate & is.null(contr)) {
    
    if (!is.null(data_sim) | !is.null(estimates_sim)) {
      
      message("Consider to provide 'contr' for your custom simulated data or analysis results.")
      
    }
    
    contr <- getContr(
      mods           = mods,
      dose_levels    = dose_levels,
      dose_weights   = n_patients,
      prior_list     = prior_list)
    
  }
  
  
  # get model names of true underlying models
  if (is.null(estimates_sim))  {
    
    model_names <- colnames(data)[-c(1, 2, 3)]
    
  } else {
    
    model_names <- "estimates_sim"
    
  }
  
  eval_design <- lapply(model_names, function (model_name) {
    
    if (is.null(estimates_sim)) {
      
      posterior_list <- getPosterior(
        data       = getModelData(data, model_name),
        prior_list = prior_list)
      
    } else {
      
      posterior_list <- mapply(getPosterior,
             mu_hat   = estimates_sim$mu_hats,
             S_hat    = estimates_sim$S_hats,
             MoreArgs = list(prior_list = prior_list),
             SIMPLIFY = FALSE)
      
      names(posterior_list) <- seq_along(estimates_sim$mu_hats)
      
    }
    
    if (reestimate & is.null(contr)) {
      
      post_sds <- sapply(posterior_list, function (post) summary(post)[, 2])
      
      contr <- apply(post_sds, 2, function (post_sd) getContr(
        mods          = mods,
        dose_levels   = dose_levels,
        cov_posterior = diag(post_sd^2)))
      
    }
    
    if (modeling) {
      
      b_mcp_mod <- performBayesianMCPMod(
        posterior_list = posterior_list,
        contr          = contr,
        crit_prob_adj  = crit_prob_adj,
        simple         = simple,
        avg_fit        = avg_fit,
        delta          = delta,
        evidence_level = evidence_level,
        med_selection  = med_selection,
        n_samples      = n_sim)
      
    } else {
      
      b_mcp <- performBayesianMCP(
        posterior_list = posterior_list,
        contr          = contr,
        crit_prob_adj  = crit_prob_adj)
      
    }
    
  })
  
  names(eval_design) <- model_names
  
  avg_success_rate <- mean(sapply(eval_design, function (x) {
    
    ifelse(identical(modeling, FALSE),
           attr(x, "successRate"),
           attr(x$BayesianMCP, "successRate"))
    
  }))
  
  attr(eval_design, "avgSuccessRate") <- avg_success_rate
  
  if (modeling & !is.null(delta)) {
    
    avg_med_id_rate <- mean(sapply(eval_design, function (x) {
      
      mean(attr(x, "MED")[, "med_reached"])
      
    }))
    
    attr(eval_design, "avgMEDIdentificationRate") <- avg_med_id_rate
    
  }
  
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