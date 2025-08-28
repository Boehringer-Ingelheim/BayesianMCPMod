#' @title getPriorList
#'
#' @param hist_data historical trial summary level data,
#' needs to be provided as a dataframe. Including information of the
#' estimates and variability.
#' @param dose_levels vector of the different doseage levels
#' @param dose_names character vector of dose levels,
#' default NULL and will be automatically created
#' based on the dose levels parameter.
#' @param robust_weight needs to be provided as a numeric
#' value for the weight of the robustification component
#'
getPriorList <- function (

  hist_data,
  dose_levels,
  dose_names       = NULL,
  robust_weight

) {

  checkmate::check_data_frame(hist_data)
  checkmate::assert_double(dose_levels, lower = 0, any.missing = FALSE)
  checkmate::check_string(dose_names, null.ok = TRUE)
  checkmate::check_vector(dose_names, null.ok = TRUE, len = length(dose_levels))
  checkmate::check_numeric(robust_weight, len = 1, null.ok = FALSE)

  sd_tot <- with(hist_data, sum(sd * n) / sum(n))

  gmap <- RBesT::gMAP(
    formula    = cbind(est, se) ~ 1 | trial,
    family     = gaussian,
    weights    = hist_data$n,
    data       = hist_data,
    beta.prior = cbind(0, 100 * sd_tot),
    tau.dist   = "HalfNormal",
    tau.prior  = cbind(0, sd_tot / 4))

  prior_ctr <- RBesT::automixfit(gmap)

  prior_ctr <- suppressMessages(RBesT::robustify(
    priormix = prior_ctr,
    weight   = robust_weight,
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

  return (prior_list)

}


# read in testdata --------------------------------------------------------

# testdata <- readRDS(file.path(getwd(), "tests/testthat/data/testdata.RDS"))
testdata <- readRDS("data/testdata.RDS")

# further setup -----------------------------------------------------------



# Create minimal test case
n_hist_trials = 2

hist_data <- data.frame(
  trial = seq(1, n_hist_trials, 1),
  est   = rep(1, n_hist_trials),
  se    = rep(1, n_hist_trials),
  sd    = rep(1, n_hist_trials),
  n     = rep(1, n_hist_trials)
)

n_patients <- c(2, 3)
dose_levels <- c(0, 2.5)
mean <- c(8, 12)
sd <- c(0.5, 0.8)

mods <- DoseFinding::Mods(
  linear = NULL,
  doses = dose_levels
)


prior_list <- getPriorList(
  hist_data   = hist_data,
  dose_levels = dose_levels,
  robust_weight = 0.5
)

n_sim = 1
alpha_crit_val = 0.05
simple = TRUE

data <- simulateData(
  n_patients  = n_patients,
  dose_levels = dose_levels,
  sd          = sd,
  mods        = mods,
  n_sim       = n_sim
)

posterior_list <- getPosterior(
  data = getModelData(data, names(mods)[1]),
  prior_list = prior_list
)

contr_mat = getContr(
  mods = mods,
  dose_levels = dose_levels,
  dose_weights = n_patients,
  prior_list = prior_list
)

crit_pval = getCritProb(
  mods = mods,
  dose_levels = dose_levels,
  dose_weights = n_patients,
  alpha_crit_val = alpha_crit_val
)

# eval_design <- assessDesign(
#   n_patients = n_patients,
#   mods = mods,
#   prior_list = prior_list,
#   n_sim = n_sim,
#   alpha_crit_val = alpha_crit_val,
#   simple = TRUE
# )

# Create covmat test case
mixnorm_test <- RBesT::mixnorm(comp1=c(0.2,2,3), comp2=c(0.2,5,6), comp3=c(0.2,8,9), comp4=c(0.2,11,12), robust=c(0.2,14,15), sigma=9.734)

mixnorm_DG1  <- RBesT::mixnorm(comp1=c(0.5,2,3), comp2=c(0.5,3,6), sigma=9.651)

mixnorm_DG2  <- RBesT::mixnorm(comp1=c(1,8,2), sigma=9.932)

mixnorm_DG3  <- RBesT::mixnorm(comp1=c(1/3,6,8), comp2=c(1/3,7,9), comp3=c(1/3,0.5,1), sigma=9.134)

mixnorm_DG4  <- RBesT::mixnorm(comp1=c(1/3,10,3), comp2=c(1/3,3,6), comp3=c(1/3,4,7), sigma=9.236)

prior_list_matrix <- vector("list", 5)
prior_list_matrix[[1]] <- mixnorm_test
prior_list_matrix[[2]] <- mixnorm_DG1
prior_list_matrix[[3]] <- mixnorm_DG2
prior_list_matrix[[4]] <- mixnorm_DG3
prior_list_matrix[[5]] <- mixnorm_DG4

names(prior_list_matrix) <- c("Ctr","DG_1","DG_2","DG_3","DG_4")

mu_hat <- c(10, 20, 30, 40, 50)
se_hat_vector <- c(1.0, 3.0, 5.0, 9.0, 6.0)
se_hat_vector_sqrt <- sqrt(se_hat_vector)

se_hat_matrix <- matrix(c(1.00, 0.00, 0.00, 0.00, 0.00,
                          0.00, 3.00, 0.00, 0.00, 0.00,
                          0.00, 0.00, 5.00, 0.00, 0.00,
                          0.00, 0.00, 0.00, 9.00, 0.00,
                          0.00, 0.00, 0.00, 0.00, 6.00), nrow = 5, ncol = 5)

se_hat_matrix2 <- matrix(c(1.00, 0.10, 0.20, 0.30, 0.40,
                           0.10, 3.00, 0.10, 0.20, 0.30,
                           0.20, 0.10, 5.00, 0.10, 0.20,
                           0.30, 0.20, 0.10, 9.00, 0.10,
                           0.40, 0.30, 0.20, 0.10, 6.00), nrow = 5, ncol = 5)

posterior <- getPosterior(
  prior_list = prior_list_matrix,
  mu_hat     = mu_hat,
  S_hat      = se_hat_matrix,
  calc_ess   = FALSE
)

posterior_noZero <- getPosterior(
  prior_list = prior_list_matrix,
  mu_hat     = mu_hat,
  S_hat      = se_hat_matrix2,
  calc_ess   = FALSE
)


# include vector based old getPosterior Function --------------------------

getPosterior_vec <- function(

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

  stopifnot("prior_list must be an object of RBesT package" =
              all(sapply(prior_list, function(x)
                methods::is(x, "normMix") |
                  methods::is(x, "betaMix") |
                  methods::is(x, "mix"))))

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

      prior_mix <- priorList2priorMix_vec(prior_list)

      posterior <- DoseFinding::mvpostmix(
        priormix = prior_mix,
        mu_hat   = mu_hat,
        S_hat    = S_hat)

      posterior_list <- postMix2posteriorList_vec(
        posterior_list = posterior,
        prior_list     = prior_list,
        calc_ess       = calc_ess)

    } else {

      posterior_list <- getPosteriorI_vec(
        prior_list = prior_list,
        mu_hat     = mu_hat,
        se_hat     = se_hat,
        calc_ess   = calc_ess)

    }

  } else if (is.null(mu_hat) && is.null(S_hat) && !is.null(data)) {

    posterior_list <- lapply(split(data, data$simulation), getPosteriorI_vec,
                             prior_list = prior_list, calc_ess = calc_ess)

  } else {

    stop ("Either 'data' or 'mu_hat' and 'se_hat' must not be NULL.")

  }

  if (length(posterior_list) == 1) {

    posterior_list <- posterior_list[[1]]

  }

  return (posterior_list)

}

getPosteriorI_vec <- function(

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

  attr(post_list, "ess") <- if (calc_ess) getESS(post_list) else numeric(0)

  diagonals <- lapply(seq_along(comp_mat_ind[, 1]), function(x) {

    lapply(seq_along(comp_mat_ind[x, ]), function(y) return(post_list[[y]]["s", comp_mat_ind[x,y]]))

  })

  attr(post_list, "full covariance matrices") <- lapply(seq_along(diagonals), function(x) diag(diagonals[[x]]))

  return (post_list)

}
priorList2priorMix_vec <- function (prior_list) {

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

postMix2posteriorList_vec <- function (

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

calcEss_vec <- function(calc_ess, posterior_output) {

  checkmate::assert_logical(calc_ess, null.ok = FALSE, len = 1)
  checkmate::assert_list(posterior_output, names = "named", any.missing = FALSE, null.ok = FALSE)

  if (calc_ess) {

    post_ESS <- getESS(posterior_output)

  } else {

    post_ESS <- numeric(0)

  }

  return(post_ESS)

}

createMapping_vec <- function (prior_list) {

  n_comps   <- unlist(lapply(prior_list, ncol))
  comp_indx <- lapply(seq_along(prior_list), function (x) seq_len(n_comps[x]))

  return (comp_indx)

}
