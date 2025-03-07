test_that("print.BayesianMCPMod works as intented", {
  expect_error(print.BayesianMCPMod())

})


test_that("print.BayesianMCP works as intented", {
  expect_error(print.BayesianMCP())

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

  b_mcp <- performBayesianMCP(
    posterior_list = posterior_list,
    contr = contr_mat,
    crit_prob_adj = crit_pval
  )

  expect_s3_class(b_mcp, "BayesianMCP")
  expect_no_error(print(b_mcp))
  expect_type(print(b_mcp), "double")

  mods <- DoseFinding::Mods(
  linear = NULL,
  emax = c(0.5, 1.2),
  exponential = 2,
  doses = c(0, 0.5, 2, 4, 8)
)
dose_levels <- c(0, 0.5, 2, 4, 8)
sd_posterior <- c(2.8, 3, 2.5, 3.5, 4)
contr_mat <- getContr(
  mods         = mods,
  dose_levels  = dose_levels,
  sd_posterior = sd_posterior
)
critVal <- getCritProb(
  mods           = mods,
  dose_weights   = c(50, 50, 50, 50, 50), # reflecting the planned sample size
  dose_levels    = dose_levels,
  alpha_crit_val = 0.6
)
prior_list <- list(
  Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 5), sigma = 2),
  DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 1, s = 12), sigma = 2),
  DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.2, s = 11), sigma = 2),
  DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.3, s = 11), sigma = 2),
  DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 2, s = 13), sigma = 2)
)
mu <- c(0, 1, 1.5, 2, 2.5)
S_hat <- c(5, 4, 6, 7, 8)
posterior_list <- getPosterior(
  prior_list = prior_list,
  mu_hat = mu,
  S_hat = S_hat
)

x <- performBayesianMCPMod(
  posterior_list = posterior_list,
  contr = contr_mat,
  crit_prob_adj = critVal,
  simple = FALSE,
  delta = 1
)

expect_s3_class(x, "BayesianMCPMod")
expect_no_error(print(x))
expect_type(print(x), "list")


})

test_that("predict.ModelFits works as intented", {
  expect_error(predict.ModelFits())

})

test_that("s3 postList functions work as intented", {
  set.seed(8080)
  dataset     <- dplyr::filter(testdata, bname == "BRINTELLIX")
  histcontrol <- dplyr::filter(dataset, dose == 0, primtime == 8, indication == "MAJOR DEPRESSIVE DISORDER",protid!=6)

  ##Create MAP Prior
  hist_data <- data.frame(
    trial = histcontrol$nctno,
    est   = histcontrol$rslt,
    se    = histcontrol$se,
    sd    = histcontrol$sd,
    n     = histcontrol$sampsize)

  dose_levels <- c(0, 2.5, 5, 10)
  post_test_list <- getPriorList(
    hist_data = hist_data,
    dose_levels = dose_levels,
    robust_weight = 0.5)
  expect_error(summary.postList())
  expect_type(summary.postList(post_test_list), "double")
  expect_error(print.postList())
  expect_type(print(post_test_list), "list")
  expect_type(print.postList(post_test_list), "list")
  # expect_true(names(print(post_test_list)) == c("Summary of Posterior Distributions",
  #                                               "Maximum Difference to Control and Dose Group",
  #                                               "Posterior Distributions"))
})

test_that("test modelFits s3 methods", {

  model_shapes <- colnames(contr_mat$contMat)
  dose_levels  <- as.numeric(rownames(contr_mat$contMat))

  model_fits  <- getModelFits(
    models      = model_shapes,
    dose_levels = dose_levels,
    posterior   = posterior_list,
    simple      = simple)

  pred <- predict(model_fits)
  pred_dosage <- predict(model_fits, doses = dose_levels)

  expect_type(pred, "list")
  expect_true(is.null(attr(pred, "doses")))
  expect_identical(attr(pred_dosage, "doses"), dose_levels)
  expect_type(print(model_fits), "list")
})
