# Test for getModelData
test_that("getModelData returns a data frame with correct columns", {

  mock_sim_data <- data.frame(
    simulation = 1:5,
    dose = rnorm(5),
    model1 = rnorm(5),
    model2 = rnorm(5)
  )
  mock_model_name <- "model1"

  result <- getModelData(sim_data = mock_sim_data, model_name = mock_model_name)

  # Assert that the result is a data frame
  expect_true(is.data.frame(result))
  # Assert that the data frame has the correct number of columns
  expect_equal(ncol(result), 3)
  # Assert that the data frame has the correct column names
  expect_equal(colnames(result), c("simulation", "dose", "response"))
})

# Test for getModelData
test_that("simulateData works as expected", {

  dataset     <- dplyr::filter(testdata, bname == "BRINTELLIX")
  histcontrol <- dplyr::filter(dataset, dose == 0, primtime == 8, indication == "MAJOR DEPRESSIVE DISORDER")

  hist_data <- data.frame(
    trial = histcontrol$nctno,
    est   = histcontrol$rslt,
    se    = histcontrol$se,
    sd    = histcontrol$sd,
    n     = histcontrol$sampsize)

  sd_tot <- with(hist_data, sum(sd * n) / sum(n))

  dose_levels <- c(0, 2.5, 5, 10, 20)

  prior_list  <- getPriorList(
    hist_data     = hist_data,
    dose_levels   = dose_levels,
    robust_weight = 0.3)


  prior_list  <- list(
    Ctr  = RBesT::mixnorm(
      comp1  = c(w = 0.446213, m = -12.774661, s = 1.393130),
      comp1  = c(w = 0.253787, m = 3.148116,   s = 3.148116),
      robust = c(w = 0.3,      m = 9.425139,   s = 9.425139),
      sigma = sd_tot),
    DG_1 = RBesT::mixnorm(
      comp1 = c(w = 1, m = -12.816875, n = 1),
      sigma = sd_tot,
      param = "mn"),
    DG_2 = RBesT::mixnorm(
      comp1 = c(w = 1, m = -12.816875, n = 1),
      sigma = sd_tot,
      param = "mn"),
    DG_3 = RBesT::mixnorm(
      comp1 = c(w = 1, m = -12.816875, n = 1),
      sigma = sd_tot,
      param = "mn"),
    DG_4 = RBesT::mixnorm(
      comp1 = c(w = 1, m = -12.816875, n = 1),
      sigma = sd_tot,
      param = "mn")
  )

  exp     <- DoseFinding::guesst(
    d     = 5,
    p     = c(0.2),
    model = "exponential",
    Maxd  = max(dose_levels))

  emax    <- DoseFinding::guesst(
    d     = 2.5,
    p     = c(0.9),
    model = "emax")

  sigemax <- DoseFinding::guesst(
    d     = c(2.5, 5),
    p     = c(0.1, 0.6),
    model = "sigEmax")

  sigemax2 <- DoseFinding::guesst(
    d     = c(2, 4),
    p     = c(0.3, 0.8),
    model = "sigEmax")

  mods <- DoseFinding::Mods(
    linear      = NULL,
    emax        = emax,
    exponential = exp,
    sigEmax     = rbind(sigemax, sigemax2),
    doses       = dose_levels,
    maxEff      = -3,
    placEff     = -12.8)

  n_patients <- c(60, 80, 80, 80, 80)

  sim_data_1 <- simulateData(
    n_patients  = n_patients,
    dose_levels = dose_levels,
    sd          = sd_tot,
    mods        = mods,
    n_sim       = 10)

  expect_type(sim_data_1, "list")

  sim_data_2 <- simulateData(
    n_patients  = n_patients,
    dose_levels = dose_levels,
    sd          = sd_tot,
    mods        = mods,
    n_sim       = 10,
    true_model = "emax")

  expect_type(sim_data_2, "list")

  sim_data_3 <- simulateData(
    n_patients  = n_patients,
    dose_levels = dose_levels,
    sd          = sd_tot,
    mods        = mods,
    n_sim       = 10,
    dr_means    = c(0.2, 0.3, 0.4, 0.5, 0.6))

  expect_type(sim_data_3, "list")

})
