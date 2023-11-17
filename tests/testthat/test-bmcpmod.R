##########################
# Tests for assessDesign #
##########################

test_that("base case input throws no error and has correct properties", {
  
  expect_no_error(
    eval_design <- assessDesign(
      n_patients = n_patients, 
      mods = mods, 
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
    attr(eval_design$linear$BayesianMCP, "dim")[1],
    n_sim
  )
  
  # assessDesign result (in this base case) should have crit_prob = 1 - alpha_crit_val
  expect_equal(
    attr(eval_design$linear$BayesianMCP, "crit_prob"),
    1 - alpha_crit_val
  )
  
})


### n_patients param ###

test_that("assessDesign validates n_patients parameter input and give appropriate error messages", {
  
  # assertions that aren't tested here for sake of brevity
    # n_patients should be a non-NULL numeric vector
  
  expect_error(
    assessDesign(n_patients = n_patients[-1], mods = mods, prior_list = prior_list, n_sim = n_sim),
    "length of n_patients should equal number of dose groups", ignore.case = T
  )
  
  expect_error(
    assessDesign(n_patients = rep(1, length(n_patients)), mods = mods, prior_list = prior_list, n_sim = n_sim),
    "at least one element in n_patients needs to be > 1", ignore.case = T
  )
})

### mods param ###

test_that("assessDesign validates mods parameter input and give appropriate error messages", {
  
  # assertions that aren't tested here for sake of brevity
    # mods should be non-NULL object of class "Mods" from {DoseFinding}
  
  
  # checking that DOseFinding didn't change how they named their 'doses' attribute
  expect_true(
    "doses" %in% names(attributes(mods))
  )
  
  mods2 <- mods
  attr(mods2, "doses") <- 0
  expect_error(
    assessDesign(n_patients = n_patients, mods = mods2, prior_list = prior_list, n_sim = n_sim),
    "number of dose groups in mods should be equal to number of dose levels in prior_list"
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
    "dose_levels" %in% names(attributes(prior_list))
  )
  
})

#########################
# Tests for getCritProb #
#########################

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

#########################
# Tests for getContrMat #
#########################

# getContrMat relies on DoseFinding, which we assumes works correctly, so the tests here are minimal

test_that("getContrMat returns the right type of object under normal case", {
  
  contr_mat = getContrMat(
    mods = mods, 
    dose_levels = dose_levels, 
    dose_weights = n_patients,
    prior_list = prior_list
  )

  expect_s3_class(
    contr_mat, "optContr"
  )

})

################################
# Tests for performBayesianMCP #
################################

test_that("performBayesianMCP returns the right type of object under normal case", {
  
  data <- simulateData(
    n_patients  = n_patients,
    dose_levels = dose_levels,
    sd          = attr(prior_list, "sd_tot"),
    mods        = mods,
    n_sim       = n_sim
  )

  posteriors_list <- getPosterior(
    data = getModelData(data, names(mods)[1]),
    prior_list = prior_list
  )

  b_mcp <- performBayesianMCP(
    posteriors_list = posteriors_list,
    contr_mat = contr_mat,
    crit_prob = crit_pval
  )

  expect_s3_class(
    b_mcp,
    "BayesianMCP"
  )

  expect_true(
    attr(b_mcp, "crit_prob") == crit_pval
  )
  
})

###################################
# Tests for performBayesianMCPMod #
###################################

test_that("performBayesianMCPMod returns the right type of object under normal case", {
  
  b_mcp_mod <- performBayesianMCPMod(
    posteriors_list = posterior_list,
    contr_mat = contr_mat,
    crit_prob = crit_pval
  )

  expect_s3_class(
    b_mcp_mod,
    "BayesianMCPMod"
  )

  expect_true(
    all(names(b_mcp_mod) == c("BayesianMCP", "Mod"))
  )
  
})


#############################
# Tests for addSignificance #
#############################

test_that("addSignificance works as intended", {
  model_fits  <- getModelFits(
    models      = colnames(contr_mat$contMat),
    dose_levels = dose_levels,
    posterior   = posteriors_list,
    simple      = TRUE
  )

  model_fits_with_sign = addSignificance(model_fits, TRUE)
  expect_true(
    model_fits_with_sign[[1]]$significant
  )

  model_fits_with_sign = addSignificance(model_fits, FALSE)
  expect_false(
    model_fits_with_sign[[1]]$significant
  )
})

#######################
# Tests for BayesMCPi #
#######################

#########################
# Tests for getPostProb #
#########################