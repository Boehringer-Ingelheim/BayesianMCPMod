#' @title getPosterior
#' 
#' @description Either the patient level data or both mu_hat as well as sd_hat must to be provided.
#' If patient level data is provided mu_hat and se_hat are calculated within the function using a linear model.
#' This function calculates the posterior for every dose group independently via the RBesT function postmix().
#' 
#' @param prior_list a prior list with information about the prior to be used for every dose group 
#' @param data dataframe containing the information of dose and response. Default NULL
#' Also a simulateData object can be provided.
#' @param mu_hat vector of estimated mean values (per dose group).
#' @param S_hat vector or matrix of estimated standard deviations (per dose group).
#' @param calc_ess boolean variable, indicating whether effective sample size should be calculated. Default FALSE
#' @return posterior_list, a posterior list object is returned with information about (mixture) posterior distribution per dose group
#' @examples
#' prior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 5), sigma = 2),
#'                    DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 1, s = 12), sigma = 2),
#'                    DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.2, s = 11), sigma = 2) ,  
#'                    DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.3, s = 11), sigma = 2) ,
#'                    DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 2, s = 13), sigma = 2))
#' mu <- c(0, 1, 1.5, 2, 2.5)
#' se <- c(5, 4, 6, 7, 8)
#' 
#' posterior_list <- getPosterior(
#'    prior_list = prior_list,
#'    mu_hat     = mu,
#'    se_hat     = se)
#'    
#' summary(posterior_list)
#' 
#' @export
getPosterior <- function(
    
  prior_list,
  data     = NULL,
  mu_hat   = NULL,
  S_hat    = NULL,
  calc_ess = FALSE
  
) {
  
  checkmate::check_data_frame(data, null.ok = TRUE)
  checkmate::check_list(prior_list, names = "named", any.missing = FALSE)
  checkmate::check_vector(mu_hat, any.missing = FALSE, null.ok = TRUE)
  checkmate::check_double(mu_hat, null.ok = TRUE, lower = -Inf, upper = Inf)
  checkmate::check_vector(S_hat, any.missing = FALSE, null.ok = TRUE)
  checkmate::check_double(S_hat, null.ok = TRUE, lower = 0, upper = Inf)
  
  if (!is.null(mu_hat) && !is.null(S_hat) && is.null(data)) {
    
    if (dim(S_hat)[2] > 1 && dim(S_hat)[2] == length(prior_list)) {
      prior_mix <- createPriorMix(prior_list)
      
      posterior <- DoseFinding::mvpostmix(
        priormix = prior_mix,
        mu_hat   = mu_hat,
        S_hat    = S_hat)
      
      posterior_list <- postmix2RBesT(
        posterior_list = posterior, 
        prior_list     = prior_list, 
        prior_mix      = prior_mix, 
        calc_ess       = calc_ess)
      
    } else if (dim(S_hat)[2] == 1) {
      
      posterior_list <- getPosteriorI(
        prior_list = prior_list,
        mu_hat     = mu_hat,
        se_hat     = S_hat,
        calc_ess   = calc_ess)
      
    }
    
  } else if (is.null(mu_hat) && is.null(S_hat) && !is.null(data)) {
    
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
  
  checkmate::check_data_frame(data_i, null.ok = TRUE)
  checkmate::check_list(prior_list, names = "named", any.missing = FALSE)
  checkmate::check_vector(mu_hat, any.missing = FALSE, null.ok = TRUE)
  checkmate::check_double(mu_hat, null.ok = TRUE, lower = -Inf, upper = Inf)
  checkmate::check_vector(se_hat, any.missing = FALSE, null.ok = TRUE)
  checkmate::check_double(se_hat, null.ok = TRUE, lower = 0, upper = Inf)
  
  if (is.null(mu_hat) && is.null(se_hat)) {
    
    checkmate::check_data_frame(data_i, null.ok = FALSE)
    checkmate::assert_names(names(data_i), must.include = "response")
    checkmate::assert_names(names(data_i), must.include = "dose")
    
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
  
  post_list <- mapply(RBesT::postmix, prior_list, m = mu_hat, se = se_hat,
                      SIMPLIFY = FALSE)
  
  if (is.null(names(prior_list))) {
    
    names(prior_list) <- c("Ctr", paste0("DG_", seq_along(post_list[-1])))
    
  }
  
  names(post_list) <- names(prior_list)
  class(post_list) <- "postList"
  
  attr(post_list, "ess") <- calcEss(calc_ess, post_list)
  
  attr(post_list, "full covariance matrices") <- replicate(length(prior_list)-1, diag(c(se_hat)), simplify = FALSE)
  #names(attr(post_list, "full covariance matrices")) <- c("comp1", "comp2", "comp3", "robust")
  
  return (post_list)
  
}

#' @title getESS
#' 
#' @description This function calculates the effective sample size for every dose group via the RBesT function ess().
#' 
#' @param post_list A posterior list object, for which the effective sample size for each dose group should be calculated
#'
#' @return A vector of the effective sample sizes for each dose group
#' 
#' @export
getESS <- function (
    
  post_list
  
) {
  
  # make s3 method for postList object
  suppressMessages(round(sapply(post_list, RBesT::ess), 1))
  
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

createPriorMix <- function(prior) {
  
  checkmate::assert_list(prior, names = "named", any.missing = FALSE, null.ok = FALSE)
  
  k <- length(prior)
  
  n_comps <- unlist(lapply(prior, ncol))
  args <- lapply(1:k, function(x) 1:n_comps[x])
  comp_ind <- do.call("expand.grid", args) 
  
  n_comps_prior <- nrow(comp_ind)
  
  prior_weight <- matrix(sapply(1:k, function(x) sapply(1:n_comps_prior, function(y) prior[[x]][1, comp_ind[y,x]])), nrow = n_comps_prior)
  prior_weight <- apply(prior_weight, 1, prod)
  prior_mean <- matrix(sapply(1:k, function(x) sapply(1:n_comps_prior, function(y) prior[[x]][2, comp_ind[y,x]])), nrow = n_comps_prior)
  prior_sd <- matrix(sapply(1:k, function(x) sapply(1:n_comps_prior, function(y) prior[[x]][3, comp_ind[y,x]])), nrow = n_comps_prior)
  
  prior_weight <- as.list(prior_weight)
  prior_mean <- asplit(prior_mean, 1)
  prior_vc <- lapply(asplit(prior_sd^2, 1), diag)
  prior_mix <- list(prior_weight, prior_mean, prior_vc)
  
  return(prior_mix)
  
}
 
postmix2RBesT <- function(
    
  posterior_list,
  prior_list,
  prior_mix,
  calc_ess = FALSE
  
) {
  
  checkmate::assert_list(posterior_list, names = "named", any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_list(prior_list, names = "named", any.missing = FALSE, null.ok = FALSE)
  
  summed_values <- vector("list", length(prior_list))
  posterior_output <- vector("list", length(prior_list))
  
  posterior_output <- lapply(seq_along(prior_list), function(i) {
    
    patterns <- lapply(prior_list[[i]]["m",], function(x) {
      
      matches <- sapply(seq_along(posterior_list$weights), function(j) {
        
        prior_mix[[2]][[j]][i] == x
        
      })
      
      result <- which(matches)
      
      return(result)
      
    })
    
    diagonals <- lapply(posterior_list$covmat, diag)
    diagonals_bound <- sapply(diagonals, cbind)
    
    summed_values <- lapply(patterns, function(x) {
      
      c(sum(as.numeric(t(t(posterior_list$weights[x])))),
        mean(as.numeric(t(sapply(posterior_list$mean, rbind))[x,i])),
        mean(diagonals_bound[i,x])
      )
      
    })
    
    additional_args <- list(sigma = sigma(prior_list[[i]]))
    
    args <- c(summed_values, additional_args)
    
    result <- do.call(RBesT::mixnorm, args)
    
    return(result)
    
  })
  
  names(posterior_output) <- c("Ctr", paste0("DG_", seq_along(posterior_output[-1])))
  
  class(posterior_output) <- "postList"
  
  attr(posterior_output, "ess") <- calcEss(calc_ess, posterior_output)
  
  attr(posterior_output, "full covariance matrices") <- posterior_list$covmat
  
  return(posterior_output)
  
}

calcEss <- function(calc_ess, posterior_output) {
  
  checkmate::assert_logical(calc_ess, null.ok = FALSE, len = 1)
  checkmate::assert_list(posterior_output, names = "named", any.missing = FALSE, null.ok = FALSE)
  
  if (calc_ess) {
    
    post_ESS <- getESS(posterior_output)
    
  } else {
    
    post_ESS <- numeric(0)
    
  }
  
  return(post_ESS)
  
}