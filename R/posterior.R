#' @title getPriorList
#' 
#' @param hist_data historical trial summary level data,
#' needs to be provided as a dataframe. Including information of the
#' estimates and variability.
#' @param dose_levels vector of the different doseage levels
#' @param dose_names character vector of dose levels,
#' default NULL and will be automatically created
#' based on the dose levels parameter.
#' @param robustify_weight Null needs to be provided as a numeric
#' value for the weight of the robustification component
#'
#' @export
getPriorList <- function (
    
  hist_data,
  dose_levels,
  dose_names       = NULL,
  robustify_weight = NULL
  
) {
  
  sd_tot <- with(hist_data, sum(sd * n) / sum(n))
  
  gmap <- RBesT::gMAP(
    formula    = cbind(est, se) ~ 1 | trial,
    weights    = hist_data$n,
    data       = hist_data,
    family     = gaussian,
    beta.prior = cbind(0, 100 * sd_tot),
    tau.dist   = "HalfNormal",
    tau.prior  = cbind(0, sd_tot / 4))
  
  prior_ctr <- RBesT::automixfit(gmap)
  
  if(is.null(robustify_weight) | !is.numeric(robustify_weight)) {
    stop("robustify_weight needs to be provided and must be numeric")
  }
  
  if (!is.null(robustify_weight)) {
    
    prior_ctr <- suppressMessages(RBesT::robustify(
      priormix = prior_ctr,
      weight   = robustify_weight,
      sigma    = sd_tot))
    
  }
  
  prior_trt <- RBesT::mixnorm(
    comp1 = c(w = 1, m = summary(prior_ctr)[1], n = 1),
    sigma = sd_tot,
    param = "mn")
  
  prior_list <- c(list(prior_ctr),
                  rep(x     = list(prior_trt),
                      times = length(dose_levels[-1])))
  
  if (is.null(dose_names)) {
    
    dose_names <- c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
    
  }
  
  names(prior_list) <- dose_names
  
  attr(prior_list, "dose_levels") <- dose_levels
  attr(prior_list, "sd_tot")      <- sd_tot
  
  return (prior_list)
  
}

#' @title getPosterior
#' 
#' @description Either the patient level data or both the mu_hat as well as the sd_hat must to be provided.
#' 
#' @param data dataframe containing the information of dose and response.
#' Also a simulateData object can be provided.
#' @param prior_list prior_list object
#' @param mu_hat vector of estimated mean values
#' @param sd_hat vector of estimated standard deviations.
#'
#' @export
getPosterior <- function(
  data,
  prior_list,
  mu_hat = NULL,
  sd_hat = NULL
  
) {
  
  posterior_list <- lapply(split(data, data$simulation), getPosteriorI,
                           prior_list = prior_list,
                           mu_hat     = mu_hat,
                           sd_hat     = sd_hat)
  
  if (length(posterior_list) == 1) {
    
    posterior_list <- posterior_list[[1]]
    
  }
  
  return (posterior_list)
  
}

getPosteriorI <- function(
    
  data_i,
  prior_list,
  mu_hat = NULL,
  sd_hat = NULL
  
) {
  
  if (is.null(mu_hat) && is.null(sd_hat)) {
    
    anova_res <- stats::lm(data_i$response ~ factor(data_i$dose) - 1)
    mu_hat    <- summary(anova_res)$coefficients[, 1]
    sd_hat    <- summary(anova_res)$coefficients[, 2]
    
  } else if (!is.null(mu_hat) && !is.null(sd_hat)) {
    
    stopifnot("m_hat length must match number of dose levels" = 
                length(prior_list) == length(mu_hat),
              "sd_hat length must match number of dose levels" = 
                length(prior_list) == length(sd_hat))
    
  } else {
    
    stop ("Both mu_hat and sd_hat must be provided.")
    
  }
  
  post_list <- mapply(RBesT::postmix, prior_list, m = mu_hat, se = sd_hat)
  
  if (is.null(names(prior_list))) {
    
    names(prior_list) <- c("Ctr", paste0("DG_", seq_along(post_list[-1])))
    
  }
  
  names(post_list) <- names(prior_list)
  class(post_list) <- "postList"
  
  return (post_list)
  
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


