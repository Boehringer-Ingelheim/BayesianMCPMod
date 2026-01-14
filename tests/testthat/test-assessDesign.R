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
  
  contr_mat <- getContr(
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
    )
  )
  
  
  sd_tot <- 9.4
  
  dose_levels <- c(0, 2.5, 5, 10, 20)
  
  prior_list <- lapply(dose_levels, function(dose_group) {
    RBesT::mixnorm(weak = c(w = 1, m = 0, s = 200), sigma = 10)
  })
  
  names(prior_list) <- c("Ctr", paste0("DG_", dose_levels[-1]))
  
  exp <- DoseFinding::guesst(
    d     = 5,
    p     = c(0.2),
    model = "exponential",
    Maxd  = max(dose_levels)
  )
  
  emax <- DoseFinding::guesst(
    d     = 2.5,
    p     = c(0.9),
    model = "emax"
  )
  
  sigemax <- DoseFinding::guesst(
    d     = c(2.5, 5),
    p     = c(0.1, 0.6),
    model = "sigEmax"
  )
  
  sigemax2 <- DoseFinding::guesst(
    d     = c(2, 4),
    p     = c(0.3, 0.8),
    model = "sigEmax"
  )
  
  mods <- DoseFinding::Mods(
    linear      = NULL,
    emax        = emax,
    exponential = exp,
    sigEmax     = rbind(sigemax, sigemax2),
    doses       = dose_levels,
    maxEff      = -3,
    placEff     = -12.8
  )
  
  n_patients <- c(60, 80, 80, 80, 80)
  
  expect_no_error(
    assessDesign(
      n_patients = n_patients,
      mods       = mods,
      prior_list = prior_list,
      sd         = sd_tot,
      n_sim      = 10,
      reestimate = TRUE
    )
  )
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