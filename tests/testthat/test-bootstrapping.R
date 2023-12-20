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
  
  result <- getBootstrapQuantiles(fit, quantiles = c(0.025,0.5, 0.975), doses = c(0, 1,2,3,4,6,8))
  expect_type(result, "list")
  
  result_2 <- getBootstrapQuantiles(fit, n_samples = 1e2, quantiles = c(0.1, 0.9), avg_fit = FALSE, doses = c(1, 2, 3))
  expect_type(result_2, "list")

  result_3 <- getBootstrapQuantiles(fit, quantiles = c(0.025,0.5, 0.975), doses = NULL)
  expect_type(result_3, "list")
})