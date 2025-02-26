# Test cases
test_that("test getBootstrapBands", {

  BMCP_result <- performBayesianMCP(
    posterior_list = posterior_list,
    contr           = contr_mat,
    crit_prob_adj   = crit_pval)

  model_shapes  <- c("linear", "emax", "exponential")

  # Option b) Making use of the complete posterior distribution
  fit <- getModelFits(
    models      = model_shapes,
    dose_levels = dose_levels,
    posterior   = posterior_list,
    simple      = FALSE)

  result <- bootstrap_quantiles <- BayesianMCPMod::getBootstrapQuantiles(
    model_fits = fit,
    quantiles  = c(0.025, 0.5, 0.975),
    doses      = c(0, 2.5, 4, 5, 7, 10),
    n_samples  = 6
  )

  expect_type(result, "list")

  result_2 <- bootstrap_quantiles <- BayesianMCPMod::getBootstrapQuantiles(
    model_fits = fit,
    quantiles  = c(0.025, 0.5, 0.975),
    doses      = c(0, 2.5, 4, 5, 7, 10),
    n_samples  = 6
  )

  expect_type(result_2, "list")
  
})
