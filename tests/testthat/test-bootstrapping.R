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
  
  bs_samples <- getBootstrapSamples(model_fits = fit,
                                    doses      = c(0, 1,2,3,4,6,8))
  
  result <- getBootstrapQuantiles(bs_samples, quantiles = c(0.025,0.5, 0.975))
  expect_type(result, "list")
  
  bs_samples_2 <- getBootstrapSamples(model_fits = fit,
                                      n_samples  = 1e2,
                                      avg_fit    = FALSE,
                                      doses      = c(1, 2, 3))
  
  result_2 <- getBootstrapQuantiles(bs_samples_2, quantiles = c(0.1, 0.9))
  expect_type(result_2, "list")
  
  bs_samples_3 <- getBootstrapSamples(model_fits = fit,
                                      doses      = NULL)

  result_3 <- getBootstrapQuantiles(bs_samples_3, quantiles = c(0.025,0.5, 0.975))
  expect_type(result_3, "list")
})