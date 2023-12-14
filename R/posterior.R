#' @title getPriorList
#' 
#' @param hist_data historical trial summary level data,
#' needs to be provided as a dataframe. Including information of the
#' estimates and variability.
#' @param dose_levels vector of the different doseage levels
#' @param dose_names character vector of dose levels,
#' default NULL and will be automatically created
#' based on the dose levels parameter.
#' @param robustify_weight needs to be provided as a numeric
#' value for the weight of the robustification component
#'
getPriorList <- function (
  
  hist_data,
  dose_levels,
  dose_names       = NULL,
  robustify_weight
  
) {
  
  checkmate::check_data_frame(hist_data)
  checkmate::assert_double(dose_levels, lower = 0, any.missing = FALSE)
  checkmate::check_string(dose_names, null.ok = TRUE)
  checkmate::check_vector(dose_names, null.ok = TRUE, len = length(dose_levels))
  checkmate::check_numeric(robustify_weight, len = 1, null.ok = FALSE)
  
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
    
    prior_ctr <- suppressMessages(RBesT::robustify(
      priormix = prior_ctr,
      weight   = robustify_weight,
      sigma    = sd_tot))

  
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
#' @param prior_list prior_list object
#' @param data dataframe containing the information of dose and response. Default NULL
#' Also a simulateData object can be provided.
#' @param mu_hat vector of estimated mean values
#' @param se_hat vector of estimated standard deviations.
#' @param calc_ess tbd. Default NULL
#'
#' @export
getPosterior <- function(
  prior_list,
  data     = NULL,
  mu_hat   = NULL,
  se_hat   = NULL,
  calc_ess = FALSE
  
) {
  checkmate::check_data_frame(data, null.ok = TRUE)
  checkmate::check_list(prior_list, names = "named", any.missing = FALSE)
  checkmate::check_vector(mu_hat, any.missing = FALSE, null.ok = TRUE)
  checkmate::check_double(mu_hat, null.ok = TRUE, lower = -Inf, upper = Inf)
  checkmate::check_vector(sd_hat, any.missing = FALSE, null.ok = TRUE)
  checkmate::check_double(sd_hat, null.ok = TRUE, lower = 0, upper = Inf)
  
  if (is.null(data)) {
    posterior_list <- getPosteriorI(data_i = NULL, prior_list = prior_list,
                                   mu_hat     = mu_hat,
                                   sd_hat     = sd_hat)
  } else {
  posterior_list <- lapply(split(data, data$simulation), getPosteriorI,
                           prior_list = prior_list,
                           mu_hat     = mu_hat,
                           sd_hat     = sd_hat)

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

  checkmate::check_data_frame(data_i, null.ok = TRUE)
  checkmate::check_list(prior_list, names = "named", any.missing = FALSE)
  checkmate::check_vector(mu_hat, any.missing = FALSE, null.ok = TRUE)
  checkmate::check_double(mu_hat, null.ok = TRUE, lower = -Inf, upper = Inf)
  checkmate::check_vector(sd_hat, any.missing = FALSE, null.ok = TRUE)
  checkmate::check_double(sd_hat, null.ok = TRUE, lower = 0, upper = Inf)
  
  if (is.null(mu_hat) && is.null(sd_hat)) {
    checkmate::check_data_frame(data_i, null.ok = FALSE)
    # checkmate::assert_names(names(data_i), must.include = "reponse")
    # needs fixing! the reponse field is not available after using simulateData
    # posterior <- getPosterior(data = simulateData(4, dose_levels, new_trial$sd, mods), prior = prior_list,
    #                           mu_hat = NULL,
    #                           sd_hat = NULL)
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
#' @description blubber
#' 
#' @param post_list blubb
#'
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


