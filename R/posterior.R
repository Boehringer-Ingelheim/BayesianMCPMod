#' @title getPosterior
#'
#' @description Either the patient level data or both mu_hat as well as S_hat must to be provided.
#' If patient level data is provided mu_hat and S_hat are calculated within the function using a linear model.
#' This function calculates the posterior distribution. Depending on the input for S_hat this step is either performed for every dose group independently via the RBesT function postmix() or the mvpostmix() function of the DoseFinding package is utilized.
#' In the latter case conjugate posterior mixture of multivariate normals are calculated (DeGroot 1970, Bernardo and Smith 1994)
#'
#' @param prior_list a prior list with information about the prior to be used for every dose group
#' @param data dataframe containing the information of dose and response. Also a simulateData object can be provided. Default NULL.
#' @param mu_hat vector of estimated mean values (per dose group). Default NULL.
#' @param S_hat covariance matrix specifying the (estimated) variability.
#' The variance-covariance matrix should be provided and the dimension of the matrix needs to match the number of dose groups. Default NULL.
#' @param calc_ess boolean variable, indicating whether effective sample size should be calculated. Default FALSE.
#' @param probability_scale A boolean to specify if the trial has a continuous or a binary outcome. Setting to TRUE will transform calculations from the logit scale to the probability scale, which can be desirable for a binary outcome. Default `attr(data, "probability_scale")`.
#' 
#' @references BERNARDO, Jl. M., and Smith, AFM (1994). Bayesian Theory. 81.
#' @return posterior_list, a posterior list object is returned with information about (mixture) posterior distribution per dose group (more detailed information about the conjugate posterior in case of covariance input for S_hat is provided in the attributes)
#' @details Kindly note that one can sample from the `posterior_list` with `lapply(posterior_list, RBesT::rmix, n = 10)`.
#' @examples
#' prior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 5), sigma = 2),
#'                    DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 1, s = 12), sigma = 2),
#'                    DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.2, s = 11), sigma = 2) ,
#'                    DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.3, s = 11), sigma = 2) ,
#'                    DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 2, s = 13), sigma = 2))
#'                    
#' mu_hat <- c(0, 1, 1.5, 2, 2.5)
#' S_hat  <- diag(c(5, 4, 6, 7, 8)^2)
#'
#' posterior_list <- getPosterior(
#'    prior_list = prior_list,
#'    mu_hat     = mu_hat,
#'    S_hat      = S_hat)
#'
#' summary(posterior_list)
#'
#' @export
getPosterior <- function(

  prior_list,
  data              = NULL,
  mu_hat            = NULL,
  S_hat             = NULL,
  calc_ess          = FALSE,
  probability_scale = attr(data, "probability_scale")

) {

  checkmate::assert_data_frame(data, null.ok = TRUE)
  checkmate::assert_list(prior_list, names = "named", any.missing = FALSE)
  checkmate::assert_vector(mu_hat, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_double(mu_hat, null.ok = TRUE, lower = -Inf, upper = Inf)
  checkmate::assert_matrix(S_hat, mode = "numeric", any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_flag(probability_scale, null.ok = TRUE)
  
  if (is.null(probability_scale)) probability_scale <- FALSE

  stopifnot("prior_list must be an object of RBesT package" =
              all(sapply(prior_list, function(x)
                methods::is(x, "normMix") |
                  methods::is(x, "betaMix") |
                  methods::is(x, "mix"))))

  if (!is.null(mu_hat) && !is.null(S_hat) && is.null(data)) {

    stopifnot("S_hat has to be a symmetrical matrix" =
              isSymmetric(S_hat))
  
    if (!isOffDiagonalZero(S_hat)) {
      
      stopifnot("dim of S_hat must equal length of prior list" =
                  dim(S_hat)[2] == length(prior_list))
      
      prior_mix <- priorList2priorMix(prior_list)
      
      post_mix <- DoseFinding::mvpostmix(
        priormix = prior_mix,
        mu_hat   = mu_hat,
        S_hat    = S_hat)
      
      posterior_list <- postMix2posteriorList(
        post_mix   = post_mix,
        prior_list = prior_list,
        calc_ess   = calc_ess)
      
    } else {
      
      ## TODO: Do we need posterior derivation based on RBesT,
      ## or should we just use DoseFinding::mvpostmix?
      ## Also interesting in terms of dependency on RBesT for validation?
      # all.equal(posterior_list2, posterior_list)
      posterior_list <- getPosteriorI(
        prior_list = prior_list,
        mu_hat     = mu_hat,
        se_hat     = diag(sqrt(S_hat)),
        calc_ess   = calc_ess)
      
    }

  } else if (is.null(mu_hat) && is.null(S_hat) && !is.null(data)) {

    posterior_list <- lapply(split(data, data$simulation), getPosteriorI,
                             prior_list        = prior_list,
                             calc_ess          = calc_ess,
                             probability_scale = probability_scale)

  } else {

    stop ("Either 'data' or 'mu_hat' and 'S_hat' must not be NULL.")

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
  calc_ess = FALSE,
  probability_scale = attr(data_i, "probability_scale")

) {

  checkmate::check_data_frame(data_i, null.ok = TRUE)
  checkmate::check_list(prior_list, names = "named", any.missing = FALSE)
  checkmate::check_vector(mu_hat, any.missing = FALSE, null.ok = TRUE)
  checkmate::check_double(mu_hat, null.ok = TRUE, lower = -Inf, upper = Inf)
  checkmate::check_vector(se_hat, any.missing = FALSE, null.ok = TRUE)
  checkmate::check_double(se_hat, null.ok = TRUE, lower = 0, upper = Inf)
  checkmate::assert_flag(probability_scale, null.ok = TRUE)
  
  if (is.null(probability_scale)) probability_scale <- FALSE

  if (is.null(mu_hat) && is.null(se_hat) && !is.null(data_i)) {

    checkmate::check_data_frame(data_i, null.ok = FALSE)
    checkmate::assert_names(names(data_i), must.include = "response")
    checkmate::assert_names(names(data_i), must.include = "dose")
    
    if (probability_scale) {
      
      logit_fit <- stats::glm(data_i$response ~ factor(data_i$dose) - 1, family = binomial)
      mu_hat    <- stats::coef(logit_fit)
      se_hat    <- diag(sqrt(stats::vcov(logit_fit)))
      
    } else {
      
      anova_res <- stats::lm(data_i$response ~ factor(data_i$dose) - 1)
      mu_hat    <- summary(anova_res)$coefficients[, 1]
      se_hat    <- summary(anova_res)$coefficients[, 2]
      
    }

  } else if (!is.null(mu_hat) && !is.null(se_hat) && is.null(data_i)) {

    stopifnot("m_hat length must match number of dose levels" =
                length(prior_list) == length(mu_hat))

  } else {

    stop ("Both mu_hat and se_hat or data_i must be provided.")

  }

  post_list <- mapply(RBesT::postmix, prior_list, m = mu_hat, se = se_hat,
                      SIMPLIFY = FALSE)

  if (is.null(names(prior_list))) {

    names(prior_list) <- c("Ctr", paste0("DG_", seq_along(post_list[-1])))

  }
  
  names(post_list) <- names(prior_list)
  class(post_list) <- "postList"

  attr(post_list, "ess") <- if (calc_ess) getESS(post_list) else numeric(0)
  
  attr(post_list, "posteriorInfo") <- priorList2priorMix(post_list)
  
  return (post_list)

}

#' @title getESS
#'
#' @description This function calculates the effective sample size for every dose group via `RBesT::ess()`.
#' @param post_list A posterior list object, for which the effective sample size for each dose group should be calculated
#' @param method A string specifying the method of ESS calculation, see `?RBesT::ess()`.
#' @param n_digits An integer for the number of digits the result should be rounded to.
#' @param ... Optional arguments applicable to specific methods, see `?RBesT::ess()`.
#' @return A vector of the effective sample sizes for each dose group
#'
#' @export
getESS <- function (

  post_list,
  method   = c("elir", "moment", "morita"),
  n_digits = 1,
  ...

) {
  
  checkmate::assert_list(post_list)
  checkmate::assert_integerish(n_digits, lower = 0, len = 1L)
  stopifnot(all(sapply(post_list, inherits, c("normMix", "mix"))))
  
  suppressMessages(round(sapply(post_list, RBesT::ess, method = method, ... = ...), 1))

}

getPostCompsI <- function (

  posterior_i

) {

  post_params <- list(
    weights = lapply(posterior_i, function (x) x[1, ]),
    means   = lapply(posterior_i, function (x) x[2, ]),
    vars    = lapply(posterior_i, function (x) x[3, ]^2))

  post_comps         <- lapply(post_params, expand.grid)
  post_comps$weights <- apply(post_comps$weights, 1, prod)

  return (post_comps)

}

priorList2priorMix <- function (prior_list) {

  checkmate::assert_list(prior_list, names = "named", any.missing = FALSE, null.ok = FALSE)

  # create mapping
  args          <- createMapping(prior_list)
  comp_ind      <- do.call("expand.grid", args)
  n_comps_prior <- nrow(comp_ind)

  # map information -> mapping function?
  prior_weight <- matrix(
    sapply(1:length(prior_list), function (x) sapply(1:n_comps_prior, function (y)
      prior_list[[x]][1, comp_ind[y, x]])), nrow = n_comps_prior)
  
  prior_mean   <- matrix(sapply(1:length(prior_list), function (x)
    sapply(1:n_comps_prior, function (y)
      prior_list[[x]][2, comp_ind[y, x]])), nrow = n_comps_prior)
  
  prior_sd     <- matrix(sapply(1:length(prior_list), function (x)
    sapply(1:n_comps_prior, function (y)
      prior_list[[x]][3, comp_ind[y, x]])), nrow = n_comps_prior)
  
  prior_weight <- apply(prior_weight, 1, prod)

  # create prior_mix object
  prior_weight <- as.list(prior_weight)
  prior_mean   <- lapply(asplit(prior_mean, 1), as.matrix)
  prior_vc     <- lapply(asplit(prior_sd^2, 1), diag)
  
  prior_mix    <- list(weights = prior_weight,
                       means   = prior_mean,
                       covMats = prior_vc)
  
  for (i in seq_len(3)) {
    
    names(prior_mix[[i]]) <- paste0("Comp", seq_along(prior_weight))
    
  }

  return (prior_mix)

}

postMix2posteriorList <- function (

  post_mix,
  prior_list,
  calc_ess = FALSE

) {

  checkmate::assert_list(post_mix, names = "named", any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_list(prior_list, names = "named", any.missing = FALSE, null.ok = FALSE)

  getIndx <- function (x, y) which(comp_mat_ind[, x] == comp_indx[[x]][y])

  # create mapping
  comp_indx    <- createMapping(prior_list)
  comp_mat_ind <- do.call("expand.grid", comp_indx)

  # map posterior information
  posterior_weight <- lapply(seq_along(prior_list), function (x)
    sapply(seq_along(comp_indx[[x]]), function (y)
      sum(unlist(post_mix$weights[getIndx(x, y)]))))

  posterior_mean <- lapply(seq_along(prior_list), function (x)
    sapply(seq_along(comp_indx[[x]]), function (y)
      mean(do.call(cbind, post_mix$mean[getIndx(x, y)])[x, ])))

  posterior_sd <- lapply(seq_along(prior_list), function (x)
    sapply(seq_along(comp_indx[[x]]), function (y)
      mean(do.call(rbind, lapply(post_mix$covmat, diag)[getIndx(x, y)])[, x])))

  combined_vectors <- mapply(function (x, y, z)
    Map(c, x, y, z), posterior_weight, posterior_mean, lapply(posterior_sd, sqrt),
    SIMPLIFY = FALSE)

  # create posterior list
  post_list <- lapply(seq_along(combined_vectors), function (x)
    do.call(RBesT::mixnorm,
            c(combined_vectors[[x]], sigma = stats::sigma(prior_list[[x]]))))

  ## fix component names
  names(post_list) <- names(prior_list)
  comp_names <- lapply(prior_list, colnames)
  
  for (i in seq_along(post_list)) {

    colnames(post_list[[i]]) <- comp_names[[i]]

  }
  rm(i)

  ## set attributes
  class(post_list) <- "postList"
  
  attr(post_list, "ess") <- if (calc_ess) getESS(post_list) else numeric(0)

  names(post_mix) <- c("weights", "means", "covMats")
  attr(post_list, "posteriorInfo") <- post_mix

  return (post_list)

}

createMapping <- function (prior_list) {

  n_comps   <- unlist(lapply(prior_list, ncol))
  comp_indx <- lapply(seq_along(prior_list), function (x) seq_len(n_comps[x]))

  return (comp_indx)

}

isOffDiagonalZero <- function (mat) {
  
  is.matrix(mat) && all(!is.na(mat)) &&
    all(abs(mat[row(mat) != col(mat)]) < .Machine$double.eps)
  
}
