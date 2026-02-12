test_that("Bootstrap samples/quantiles: probability_scale=TRUE returns bounded values and stable quantile handling", {
  skip_if_not_installed("RBesT")
  skip_if_not_installed("DoseFinding")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  
  set.seed(12)
  
  dose_levels <- c(0, 1, 2, 4)
  probs <- c(0.10, 0.20, 0.30, 0.55)
  
  # Build a posterior directly from mu_hat/se_hat to keep test fast and stable
  mu_hat <- qlogis(probs)
  se_hat <- rep(0.35, length(dose_levels))
  
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = mu_hat[i], s = 1.5), sigma = 2)
    }),
    c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
  )
  
  post <- getPosterior(
    prior_list        = prior_list,
    mu_hat            = mu_hat,
    S_hat             = diag(se_hat^2),
    probability_scale = TRUE
  )
  
  mods <- DoseFinding::Mods(
    linear = NULL,
    emax   = c(0.5, 1.2),
    doses  = dose_levels,
    maxEff = 2
  )
  
  fits <- getModelFits(
    models            = mods,
    dose_levels       = dose_levels,
    posterior         = post,
    simple            = TRUE,
    avg_fit           = TRUE,
    probability_scale = TRUE
  )
  
  bs <- getBootstrapSamples(
    model_fits        = fits,
    n_samples         = 60,
    doses             = dose_levels,
    probability_scale = TRUE
  )
  
  expect_true(all(c("model", "dose", "sample_id", "abs", "diff") %in% names(bs)))
  expect_true(all(bs$abs >= 0 & bs$abs <= 1))
  # diff is probability difference vs control => in [-1,1]
  expect_true(all(bs$diff >= -1 & bs$diff <= 1))
  
  # Quantiles: duplicates + unsorted are handled via sort(unique())
  qs_in <- c(0.5, 0.025, 0.5, 0.975)
  bq <- getBootstrapQuantiles(
    model_fits        = fits,
    quantiles         = qs_in,
    n_samples         = 60,
    doses             = dose_levels,
    probability_scale = TRUE
  )
  
  expect_true(all(c("model", "dose", "q_prob", "q_val", "sample_type") %in% names(bq)))
  expect_equal(sort(unique(qs_in)), sort(unique(bq$q_prob)))
  expect_true(all(bq$sample_type %in% c("abs", "diff")))
  
  # Bounds by sample_type
  expect_true(all(bq$q_val[bq$sample_type == "abs"]  >= 0 & bq$q_val[bq$sample_type == "abs"] <= 1))
  expect_true(all(bq$q_val[bq$sample_type == "diff"] >= -1 & bq$q_val[bq$sample_type == "diff"] <= 1))
})
