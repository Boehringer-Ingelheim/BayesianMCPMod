# Test predictModelFit function
test_that("predictModelFit works correctly", {
  model_fit <- list(
    model = "emax",
    coeffs = c(e0 = 0, eMax = 1, ed50 = 2),
    dose_levels = c(0, 1, 2, 3, 4, 5)
  )

  pred_vals <- predictModelFit(model_fit)
  expect_type(pred_vals, "double")

  dose_levels <- c(0, 2.5, 5, 10)

  prior_list <- lapply(dose_levels, function(dose_group) {
    RBesT::mixnorm(weak = c(w = 1, m = 0, s = 200), sigma = 10)
  })

  names(prior_list) <- c("Ctr", paste0("DG_", dose_levels[-1]))

  trial_data <- dplyr::filter(
    dplyr::filter(tibble::tibble(testdata), bname == "BRINTELLIX"),
    primtime == 8,
    indication == "MAJOR DEPRESSIVE DISORDER",
    protid == 5
  )

  posterior <- getPosterior(
    prior_list = prior_list,
    mu_hat     = trial_data$rslt,
    S_hat      = trial_data$se,
    calc_ess = TRUE
  )

  n_patients <- c(128, 124, 129, 122)

  # Guesstimate estimation
  exp_guesst <- DoseFinding::guesst(
    model = "exponential",
    d = 5, p = 0.2, Maxd = max(dose_levels)
  )
  emax_guesst <- DoseFinding::guesst(
    model = "emax",
    d = 2.5, p = 0.9
  )
  sigEmax_guesst <- DoseFinding::guesst(
    model = "sigEmax",
    d = c(2.5, 5), p = c(0.5, 0.95)
  )
  logistic_guesst <- DoseFinding::guesst(
    model = "logistic",
    d = c(5, 10), p = c(0.1, 0.85)
  )

  betaMod_params <- c(delta1 = 1, delta2 = 1)
  quadratic_params <- c(delta2 = -0.1)

  mods <- DoseFinding::Mods(
    linear      = NULL,
    # guesstimate scale
    exponential = exp_guesst,
    emax        = emax_guesst,
    sigEmax     = sigEmax_guesst,
    logistic    = logistic_guesst,
    # parameter scale
    betaMod     = betaMod_params,
    quadratic   = quadratic_params,
    # Options for all models
    doses       = dose_levels,
    maxEff      = -1,
    placEff     = -12.8
  )

  fit_simple <- getModelFits(
    models      = names(mods),
    dose_levels = dose_levels,
    posterior   = posterior,
    simple      = TRUE
  )

  expect_no_error(fit_simple)
  expect_type(fit_simple, "list")

  fit <- getModelFits(
    models      = names(mods),
    dose_levels = dose_levels,
    posterior   = posterior,
    simple      = FALSE
  )

  expect_no_error(fit)
  expect_type(fit, "list")
})

# Test for addModelWeights
test_that("addModelWeights works correctly in a simple case", {
  # Create a test case
  test_model_fits1 <- list(
    model1 = list(gAIC = 1),
    model2 = list(gAIC = 1),
    model3 = list(gAIC = 1)
  )
  expected_output1 <- list(
    list(gAIC = 1, model_weight = 1 / 3),
    list(gAIC = 1, model_weight = 1 / 3),
    list(gAIC = 1, model_weight = 1 / 3)
  )

  # Call the function with the test case
  result <- addModelWeights(model_fits = test_model_fits1)
  # Assert the expected behavior
  expect_true(is.list(result))
  expect_equal(length(result), length(test_model_fits1))
  expect_true(all(sapply(result, function(x) "model_weight" %in% names(x))))
  expect_equal(result, expected_output1)

  test_model_fits2 <- list(
    model1 = list(gAIC = 100),
    model2 = list(gAIC = 0),
    model3 = list(gAIC = 100)
  )
  expected_output2 <- list(
    list(gAIC = 100, model_weight = 0),
    list(gAIC = 0, model_weight = 1),
    list(gAIC = 100, model_weight = 0)
  )
  result <- addModelWeights(model_fits = test_model_fits2)
  expect_equal(result, expected_output2, tolerance = 1e-10)
})

# Test for getGenAIC
test_that("getGenAIC calculates AIC correctly using snapshot test and a simple example", {
  # Create a test case
  test_model_fit <- list(
    pred_values = rep(1, 3),
    coeffs = rep(1, 3)
  )
  test_post_combs <- list(
    means = matrix(rep(1, 6), nrow = 2),
    vars = matrix(rep(1, 6), nrow = 2),
    weights = c(1, 1)
  )

  result <- getGenAIC(model_fit = test_model_fit, post_combs = test_post_combs)
  # expected result determined with snapshot of behavior prior to first release on CRAN
  expected_result <- 6
  # Assert the expected behavior
  expect_equal(result, expected_result)
})
