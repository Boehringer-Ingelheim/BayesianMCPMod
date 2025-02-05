
# Tests for assessDesign --------------------------------------------------



test_that("base case input throws no error and has correct properties", {

  expect_no_error(
    eval_design <- assessDesign(
      n_patients = n_patients,
      mods = mods,
      sd = sd,
      prior_list = prior_list,
      n_sim = n_sim,
      alpha_crit_val = alpha_crit_val,
      simple = TRUE
    )
  )

  # assessDesign should give results for each model in mods
  expect_equal(
    names(eval_design), names(mods)
  )

  # assessDesign result should have rows = n_sim
  expect_equal(
    attr(eval_design$linear, "dim")[1],
    n_sim
  )

  # assessDesign result (in this base case) should have crit_prob = 1 - alpha_crit_val
  expect_equal(
    attr(eval_design$linear, "critProb"),
    1 - alpha_crit_val
  )

  contr_mat = getContr(
    mods = mods,
    dose_levels = dose_levels,
    dose_weights = n_patients,
    prior_list = prior_list
  )
  
  expect_no_error(
    eval_design <- assessDesign(
      n_patients = n_patients,
      mods = mods,
      sd = sd,
      prior_list = prior_list,
      n_sim = n_sim,
      alpha_crit_val = alpha_crit_val,
      simple = TRUE,
      modeling = TRUE
    )
  )
  
  # assessDesign result should have rows = n_sim
  expect_equal(
    attr(eval_design$linear$BayesianMCP, "dim")[1],
    n_sim
  )
  
  # assessDesign result (in this base case) should have crit_prob = 1 - alpha_crit_val
  expect_equal(
    attr(eval_design$linear$BayesianMCP, "critProb"),
    1 - alpha_crit_val
  )

  expect_no_error(
    assessDesign(
    n_patients = n_patients,
    mods = mods,
    sd = sd,
    prior_list = prior_list,
    n_sim = n_sim,
    alpha_crit_val = alpha_crit_val,
    simple = TRUE,
    reestimate = TRUE,
    contr = contr_mat
  ))


  sd_tot <- 9.4

  dose_levels <- c(0, 2.5, 5, 10, 20)

  prior_list <- lapply(dose_levels, function(dose_group) {
    RBesT::mixnorm(weak = c(w = 1, m = 0, s = 200), sigma = 10)
  })

  names(prior_list) <- c("Ctr", paste0("DG_", dose_levels[-1]))

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

  expect_no_error(
    assessDesign(
      n_patients  = n_patients,
      mods        = mods,
      prior_list  = prior_list,
      sd          = sd_tot,
      n_sim       = 10,
      reestimate = TRUE
      ))

})


### n_patients param ###

test_that("assessDesign validates n_patients parameter input and give appropriate error messages", {

  # assertions that aren't tested here for sake of brevity
    # n_patients should be a non-NULL numeric vector

  expect_error(
    assessDesign(n_patients = n_patients[-1], sd = sd, mods = mods, prior_list = prior_list, n_sim = n_sim)
  )

  expect_error(
    assessDesign(n_patients = rep(1, length(n_patients)), sd = sd, mods = mods, prior_list = prior_list, n_sim = n_sim),
  )
  
})

### mods param ###

test_that("assessDesign validates mods parameter input and give appropriate error messages", {

  # assertions that aren't tested here for sake of brevity
    # mods should be non-NULL object of class "Mods" from {DoseFinding}


  # checking that DoseFinding didn't change how they named their 'doses' attribute
  expect_true(
    "doses" %in% names(attributes(mods))
  )

  mods2 <- mods
  attr(mods2, "doses") <- 0
  expect_error(
    assessDesign(n_patients = n_patients, mods = mods2, sd = sd, prior_list = prior_list, n_sim = n_sim)
  )
  rm(mods2)

})

## prior_list param ###

test_that("assessDesign validates prior_list parameter input and give appropriate error messages", {

  # assertions that aren't tested here for sake of brevity
    # prior_list should be a non-NULL named list with length = number of dose levels
    # length(attr(prior_list, "dose_levels")) == n_patients (see above)

  # checking that we didn't change how we named the 'dose_levels' attribute
  expect_true(
    "doses" %in% names(attributes(mods))
  )

})


# Tests for getCritProb ---------------------------------------------------



# getCritProb relies on DoseFinding, which we assumes works correctly, so the tests here are minimal

test_that("getCritProb returns the right type of value under normal case", {

  crit_pval = getCritProb(
    mods = mods,
    dose_levels = dose_levels,
    dose_weights = n_patients,
    alpha_crit_val = alpha_crit_val
  )

  expect_type(
    crit_pval, "double"
  )

  expect_true(
    crit_pval >= 0 & crit_pval <= 1
  )

})


# Tests for getContrMat ---------------------------------------------------



# getContrMat relies on DoseFinding, which we assumes works correctly, so the tests here are minimal

test_that("getContrMat returns the right type of object under normal case", {

  contr_mat = getContr(
    mods = mods,
    dose_levels = dose_levels,
    dose_weights = n_patients,
    prior_list = prior_list
  )

  expect_s3_class(
    contr_mat, "optContr"
  )

})

test_that("getContrMat works as expected", {

  dose_levels <- c(0, 2.5, 5, 10)
  sd_posterior <- c(2.8, 3, 2.5, 3.5)

  contr_mat_post_sd <- getContr(
    mods         = mods,
    dose_levels  = dose_levels,
    sd_posterior = sd_posterior
  )

  se_new_trial <- c(0.3, 0.7, 0.9, 2.1)

  contr_mat_se_new = getContr(
    mods = mods,
    dose_levels = dose_levels,
    se_new_trial = se_new_trial
  )


  expect_s3_class(
    contr_mat_post_sd, "optContr"
  )

  expect_no_error(contr_mat_post_sd)

  expect_s3_class(
    contr_mat_se_new, "optContr"
  )

  expect_no_error(contr_mat_se_new)


  expect_error(
   getContr(
      mods = mods,
      dose_levels = dose_levels
    )
  )

})


# Tests for performBayesianMCP --------------------------------------------



test_that("performBayesianMCP returns the right type of object under normal case", {

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

  expect_s3_class(
    b_mcp,
    "BayesianMCP"
  )

  expect_true(
    attr(b_mcp, "critProbAdj") == crit_pval
  )

  expect_type(
    attr(b_mcp, "essAvg"), "logical"
  )

  expect_type(
    attr(b_mcp, "successRate"), "double"
  )


  expect_type(b_mcp, "double")

})


# Tests for performBayesianMCPMod -----------------------------------------



test_that("performBayesianMCPMod returns the right type of object under normal case", {

  b_mcp_mod <- performBayesianMCPMod(
    posterior_list = posterior_list,
    contr = contr_mat,
    crit_prob_adj = crit_pval
  )

  expect_s3_class(
    b_mcp_mod,
    "BayesianMCPMod"
  )

  expect_true(
    all(names(b_mcp_mod) == c("BayesianMCP", "Mod"))
  )

})


# Tests for addSignificance -----------------------------------------------



test_that("addSignificance works as intended", {
  model_fits  <- list(linear = 1)

  model_fits_with_sign = addSignificance(model_fits, list(TRUE))
  expect_true(
    model_fits_with_sign[[1]][["significant"]]
  )

  model_fits_with_sign = addSignificance(model_fits, list(FALSE))
  expect_false(
    model_fits_with_sign[[1]]$significant
  )
})



# Tests for getPostProb ---------------------------------------------------

# Test for getPostProb
test_that("getPostProb works correctly in a simple case", {
  # Create a test case
  contr_j <- c(1, 1)
  post_combs_i <- list(
    means = matrix(c(0, 0, 1, 1), nrow = 2),
    vars = matrix(c(0, 0, 1, 1), nrow = 2),
    weights = c(0.5, 0.5)
  )
  # Call the function with the test case
  result <- getPostProb(
    contr_j,
    post_combs_i
  )
  # Assert the expected behavior
  expect_equal(result, stats::pnorm(1))
})
