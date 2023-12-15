
#' @title getPosterior
#' 
#' @description Either the patient level data or both mu_hat as well as sd_hat must to be provided. If patient level data is provided mu_hat and se_hat are calculated within the function using a linear model.
#' This function calculates the posterior for every dose group independently via the RBesT function postmix.
#' 
#' @param prior_list prior_list object
#' @param data dataframe containing the information of dose and response. Default NULL
#' Also a simulateData object can be provided.
#' @param mu_hat vector of estimated mean values (per dose group).
#' @param se_hat vector of estimated standard deviations (per dose group).
#' @param calc_ess boolean variable, indicating whether effective sample size should be calculated. Default FALSE
#' @return posterior_list, a posterior list object is returned with information about (mixture) posterior distribution per dose group
#' @export
getPosterior <- function(
  prior_list,
  data     = NULL,
  mu_hat   = NULL,
  se_hat   = NULL,
  calc_ess = FALSE
  
) {
  
  if (!is.null(mu_hat) && !is.null(se_hat) && is.null(data)) {
    
    posterior_list <- getPosteriorI(
      prior_list = prior_list,
      mu_hat     = mu_hat,
      se_hat     = se_hat,
      calc_ess   = calc_ess)
    
  } else if (is.null(mu_hat) && is.null(se_hat) && !is.null(data)) {
    
    posterior_list <- lapply(split(data, data$simulation), getPosteriorI,
                             prior_list = prior_list, calc_ess = calc_ess)
    
  } else {
    
    stop ("Either 'data' or 'mu_hat' and 'se_hat' must not be NULL.")
    
  }
 
  if (length(posterior_list) == 1) {
    
    posterior_list <- posterior_list[[1]]
    
  }
  
  return (posterior_list)
  
}

getPosteriorI <- function(
    
  data_i   = NULL,
  prior_list,
  mu_hat   = NULL,
  se_hat   = NULL,
  calc_ess = FALSE
  
) {
  
  if (is.null(mu_hat) && is.null(se_hat)) {
    
    anova_res <- stats::lm(data_i$response ~ factor(data_i$dose) - 1)
    mu_hat    <- summary(anova_res)$coefficients[, 1]
    se_hat    <- summary(anova_res)$coefficients[, 2]
    
  } else if (!is.null(mu_hat) && !is.null(se_hat)) {
    
    stopifnot("m_hat length must match number of dose levels" = 
                length(prior_list) == length(mu_hat),
              "se_hat length must match number of dose levels" = 
                length(prior_list) == length(se_hat))
    
  } else {
    
    stop ("Both mu_hat and se_hat must be provided.")
    
  }
  
  post_list <- mapply(RBesT::postmix, prior_list, m = mu_hat, se = se_hat)
  
  if (is.null(names(prior_list))) {
    
    names(prior_list) <- c("Ctr", paste0("DG_", seq_along(post_list[-1])))
    
  }
  
  names(post_list)       <- names(prior_list)
  class(post_list)       <- "postList"
  attr(post_list, "ess") <- ifelse(
    test = calc_ess,
    yes  = getESS(post_list),
    no   = numeric(0))
  
  return (post_list)
  
}

#' @title getESS
#' 
#' @description This function calculates the effective sample size for every dose group via the RBesT function ess.
#' 
#' @param post_list a posterior list object, for which the effective sample size (per dose group) should be calculated
#'
#' @return a vector of the effective sample sizes (per dose group) 
#' @export
getESS <- function (
    
  post_list
  
) {
  
  suppressMessages(sapply(post_list, RBesT::ess))
  
}

getPostCombsI <- function (
    
  posterior_i
  
) {
  
  post_params <- list(
    weights = lapply(posterior_i, function (x) x[1, ]),
    means   = lapply(posterior_i, function (x) x[2, ]),
    vars    = lapply(posterior_i, function (x) x[3, ]^2))
  
  post_combs         <- lapply(post_params, expand.grid)
  post_combs$weights <- apply(post_combs$weights, 1, prod)
  
  return (post_combs)
  
}


