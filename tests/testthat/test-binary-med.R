test_that("getMED: probability scale enforces delta bounds and returns expected structure", {
  skip_if_not_installed("RBesT")
  skip_if_not_installed("DoseFinding")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  
  set.seed(13)
  
  dose_levels <- c(0, 1, 2, 4)
  probs <- c(0.10, 0.18, 0.28, 0.45)
  
  mu_hat <- qlogis(probs)
  se_hat <- rep(0.40, length(dose_levels))
  
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = mu_hat[i], s = 1.2), sigma = 2)
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
  
  # delta bounds on probability scale
  expect_error(getMED(delta = 1.5, model_fits = fits, probability_scale = TRUE))
  expect_error(getMED(delta = -0.1, model_fits = fits, probability_scale = TRUE))
  
  med <- getMED(delta = 0.10, model_fits = fits, probability_scale = TRUE)
  expect_true(is.matrix(med))
  expect_equal(rownames(med), c("med_reached", "med"))
  expect_true(all(colnames(med) %in% names(fits)))
  
  # Bootstrap path
  bq <- getBootstrapQuantiles(
    model_fits        = fits,
    quantiles         = c(0.2, 0.8),
    n_samples         = 80,
    doses             = dose_levels,
    probability_scale = TRUE
  )
  
  med2 <- getMED(
    delta             = 0.10,
    evidence_level    = 0.8,
    bs_quantiles      = bq,
    probability_scale = TRUE
  )
  expect_true(is.matrix(med2))
  expect_equal(rownames(med2), c("med_reached", "med"))
})
