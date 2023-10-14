#' @title getPosterior
#' 
#' @param data tbd
#' @param prior_list prior_list
#' @param mu_hat tbd
#' @param sd_hat tbd
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
  
  names(post_list) <- c("Ctr", paste0("DG_", seq_along(post_list[-1])))
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


