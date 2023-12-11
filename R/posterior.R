#' @title getPosterior
#' 
#' @param data tbd
#' @param prior_list prior_list
#' @param mu_hat tbd
#' @param se_hat tbd
#'
#' @export
getPosterior <- function(
    
  prior_list,
  data   = NULL,
  mu_hat = NULL,
  se_hat = NULL
  
) {
  
  if (!is.null(mu_hat) && !is.null(se_hat) && is.null(data)) {
    
    posterior_list <- getPosteriorI(
      prior_list = prior_list,
      mu_hat     = mu_hat,
      se_hat     = se_hat)
    
  } else if (is.null(mu_hat) && is.null(se_hat) && !is.null(data)) {
    
    posterior_list <- lapply(split(data, data$simulation), getPosteriorI,
                             prior_list = prior_list)
    
  } else {
    
    stop ("Either 'data' or 'mu_hat' and 'se_hat' must not be NULL.")
    
  }
  
  if (length(posterior_list) == 1) {
    
    posterior_list <- posterior_list[[1]]
    
  }
  
  return (posterior_list)
  
}

getPosteriorI <- function(
    
  data_i = NULL,
  prior_list,
  mu_hat = NULL,
  se_hat = NULL
  
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


