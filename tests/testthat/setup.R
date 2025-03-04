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

testdata <- readRDS("data/testdata.RDS")



# further setup -----------------------------------------------------------



getPostProb <- function (

  contr_j,     # j: dose level
  post_combs_i # i: simulation outcome

) {

  ## Test statistic = sum over all components of
  ## posterior weight * normal probability distribution of
  ## critical values for doses * estimated mean / sqrt(product of critical values for doses)

  ## Calculation for each component of the posterior
  contr_theta   <- apply(post_combs_i$means, 1, `%*%`, contr_j)
  contr_var     <- apply(post_combs_i$vars, 1, `%*%`, contr_j^2)
  contr_weights <- post_combs_i$weights

  ## P(c_m * theta > 0 | Y = y) for a shape m (and dose j)
  post_probs <- sum(contr_weights * stats::pnorm(contr_theta / sqrt(contr_var)))

  return (post_probs)

}

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
se_hat_vector_sqrt <- c(sqrt(1), sqrt(3), sqrt(5), sqrt(9), sqrt(6))

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
