#' @title getPosterior
#'
#' @description Either the patient level data or both mu_hat as well as S_hat must to be provided.
#' If patient level data is provided mu_hat and S_hat are calculated within the function using a linear model.
#' This function calculates the posterior distribution. Depending on the input for S_hat this step is either performed for every dose group independently via the RBesT function postmix() or the mvpostmix() function of the DoseFinding package is utilized.
#' In the latter case conjugate posterior mixture of multivariate normals are calculated (DeGroot 1970, Bernardo and Smith 1994)
#'
#' @param prior_list a prior list with information about the prior to be used for every dose group
#' @param data dataframe containing the information of dose and response. Default NULL
#' Also a simulateData object can be provided.
#' @param mu_hat vector of estimated mean values (per dose group).
#' @param S_hat Either a vector or a covariance matrix specifying the (estimated) variability can be specified. The length of the vector (resp. the dimension of the matrix) needs to match the number of dose groups. Please note that for a vector input the numbers should reflect the standard error per dose group (i.e. square root of variance), while for a matrix input the variance-covariance matrix should be provided.
#' @param calc_ess boolean variable, indicating whether effective sample size should be calculated. Default FALSE
#' @references BERNARDO, Jl. M., and Smith, AFM (1994). Bayesian Theory. 81.
#' @return posterior_list, a posterior list object is returned with information about (mixture) posterior distribution per dose group (more detailed information about the conjugate posterior in case of covariance input for S_hat is provided in the attributes)
#' @examples
#' prior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 5), sigma = 2),
#'                    DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 1, s = 12), sigma = 2),
#'                    DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.2, s = 11), sigma = 2) ,
#'                    DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.3, s = 11), sigma = 2) ,
#'                    DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 2, s = 13), sigma = 2))
#' mu <- c(0, 1, 1.5, 2, 2.5)
#' S_hat <- c(5, 4, 6, 7, 8)
#'
#' posterior_list <- getPosterior(
#'    prior_list = prior_list,
#'    mu_hat     = mu,
#'    S_hat     = S_hat)
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
  
  is_matrix_S_hat <- FALSE
  
  stopifnot("prior_list must be an object of RBesT package" =
              all(sapply(prior_list, function(x) methods::is(x, "normMix") |
                           methods::is(x, "betaMix") | methods::is(x, "mix"))))

  if (!is.null(mu_hat) && !is.null(S_hat) && is.null(data)) {
    
    if (is.matrix(S_hat)) {
      
      is_matrix_S_hat <- TRUE
      
    } else if (is.vector(S_hat)) {
      
      is_matrix_S_hat <- FALSE
      
      se_hat <- S_hat
      rm(S_hat)
      
    } else {
      
      stop("S_hat has to be either a vector or matrix")
      
    }
    
    if (is_matrix_S_hat) {
      
      stopifnot("dim of S_hat must equal length of prior list" =
                  dim(S_hat)[2] == length(prior_list))
      
      prior_mix <- priorList2priorMix(prior_list)
      
      posterior <- DoseFinding::mvpostmix(
        priormix = prior_mix,
        mu_hat   = mu_hat,
        S_hat    = S_hat)
      
      posterior_list <- postMix2posteriorList(
        posterior_list = posterior,
        prior_list     = prior_list,
        calc_ess       = calc_ess)
      
    } else if (!is_matrix_S_hat) {
      
      posterior_list <- getPosteriorI(
        prior_list = prior_list,
        mu_hat     = mu_hat,
        se_hat     = se_hat,
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
                length(prior_list) == length(mu_hat))
    # ,
    #           "se_hat length must match number of dose levels" =
    #             length(prior_list) == length(se_hat))
    
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
  
  comp_indx    <- createMapping(prior_list)
  comp_mat_ind <- do.call("expand.grid", comp_indx)
  
  attr(post_list, "ess") <- calcEss(calc_ess, post_list)
  
  diagonals <- lapply(seq_along(comp_mat_ind[, 1]), function(x) {
    
    lapply(seq_along(comp_mat_ind[x, ]), function(y) return(post_list[[y]]["s", comp_mat_ind[x,y]]))
    
  })
  
  attr(post_list, "full covariance matrices") <- lapply(seq_along(diagonals), function(x) diag(diagonals[[x]]))
  
  return (post_list)
  
}

#' @title getESS
#'
#' @description This function calculates the effective sample size for every dose group via the RBesT function ess().
#' @param post_list A posterior list object, for which the effective sample size for each dose group should be calculated
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

priorList2priorMix <- function (prior_list) {
  
  checkmate::assert_list(prior_list, names = "named", any.missing = FALSE, null.ok = FALSE)
  
  # create mapping
  args <- createMapping(prior_list)
  comp_ind <- do.call("expand.grid", args)
  n_comps_prior <- nrow(comp_ind)
  
  # map information -> mapping function?
  prior_weight <- matrix(
    sapply(1:length(prior_list), function (x) sapply(1:n_comps_prior,
                                                     function (y) prior_list[[x]][1, comp_ind[y, x]])), nrow = n_comps_prior)
  
  prior_mean   <- matrix(sapply(1:length(prior_list), function (x) sapply(1:n_comps_prior, function (y) prior_list[[x]][2, comp_ind[y, x]])), nrow = n_comps_prior)
  prior_sd     <- matrix(sapply(1:length(prior_list), function (x) sapply(1:n_comps_prior, function (y) prior_list[[x]][3, comp_ind[y, x]])), nrow = n_comps_prior)
  
  prior_weight <- apply(prior_weight, 1, prod)
  
  # create prior_mix object
  prior_weight <- as.list(prior_weight)
  prior_mean   <- asplit(prior_mean, 1)
  prior_vc     <- lapply(asplit(prior_sd^2, 1), diag)
  prior_mix    <- list(prior_weight, prior_mean, prior_vc)
  
  return (prior_mix)
  
}

postMix2posteriorList <- function (
    
  posterior_list,
  prior_list,
  calc_ess = FALSE
  
) {
  
  checkmate::assert_list(posterior_list, names = "named", any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_list(prior_list, names = "named", any.missing = FALSE, null.ok = FALSE)
  
  getIndx       <- function (x, y) which(comp_mat_ind[, x] == comp_indx[[x]][y])
  
  # create mapping
  comp_indx    <- createMapping(prior_list)
  comp_mat_ind <- do.call("expand.grid", comp_indx)
  
  # map posterior information
  posterior_weight <- lapply(seq_along(prior_list), function (x)
    sapply(seq_along(comp_indx[[x]]), function (y)
      sum(unlist(posterior_list$weights[getIndx(x, y)]))))
  
  posterior_mean <- lapply(seq_along(prior_list), function (x)
    sapply(seq_along(comp_indx[[x]]), function (y)
      mean(do.call(cbind, posterior_list$mean[getIndx(x, y)])[x, ])))
  
  posterior_sd <- lapply(seq_along(prior_list), function (x)
    sapply(seq_along(comp_indx[[x]]), function (y)
      mean(do.call(rbind, lapply(posterior_list$covmat, diag)[getIndx(x, y)])[, x])))
  
  combined_vectors <- mapply(function (x, y, z)
    Map(c, x, y, z), posterior_weight, posterior_mean, lapply(posterior_sd, sqrt),
    SIMPLIFY = FALSE)
  
  # create posterior list
  posterior_list_RBesT <- lapply(seq_along(combined_vectors), function (x)
    do.call(RBesT::mixnorm,
            c(combined_vectors[[x]], sigma = stats::sigma(prior_list[[x]]))))

  ## fix component names
  names(posterior_list_RBesT) <- names(prior_list)
  comp_names <- lapply(prior_list, colnames)
  for (i in seq_along(posterior_list_RBesT)) {
    
    colnames(posterior_list_RBesT[[i]]) <- comp_names[[i]]
    
  }
  rm(i)
  
  ## set attributes
  class(posterior_list_RBesT) <- "postList"
  attr(posterior_list_RBesT, "ess") <- calcEss(calc_ess, posterior_list_RBesT)
  
  diagonals <- lapply(seq_along(comp_mat_ind[, 1]), function(x) {
    
    lapply(seq_along(comp_mat_ind[x, ]), function(y) return(posterior_list_RBesT[[y]]["s", comp_mat_ind[x,y]]))
    
  })
  
  attr(posterior_list_RBesT, "covariance matrices") <- lapply(seq_along(diagonals), function(x) diag(diagonals[[x]]))
  
  return (posterior_list_RBesT)
  
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

createMapping <- function (prior_list) {
  
  n_comps   <- unlist(lapply(prior_list, ncol))
  comp_indx <- lapply(seq_along(prior_list), function (x) seq_len(n_comps[x]))
  
  return (comp_indx)
  
}
