##########################
# Tests for assessDesign #
##########################

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
    assessDesign(n_patients = n_patients[-1], sd = sd, mods = mods, prior_list = prior_list, n_sim = n_sim),
    "length of n_patients should equal number of dose groups", ignore.case = T
  )
  
  expect_error(
    assessDesign(n_patients = rep(1, length(n_patients)), sd = sd, mods = mods, prior_list = prior_list, n_sim = n_sim),
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
    assessDesign(n_patients = n_patients, mods = mods2, sd = sd, prior_list = prior_list, n_sim = n_sim),
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

  posterior_list <- getPosterior(
    data = getModelData(data, names(mods)[1]),
    prior_list = prior_list
  )

  b_mcp <- performBayesianMCP(
    posterior_list = posterior_list,
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
    posterior_list = posterior_list,
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

#######################
# Tests for BayesMCPi #
#######################

test_that("BayesMCPi function works correctly in a simple case", {

  # BayesMCPi should return a list of length 3 named with "sign", "p_val", and "post_probs"
  # The logic being tested here is: 
    # BayesMCPi returns 1 if the posterior probability is strictly greater than the critical value, and 0 otherwise
  
  # Define a mock function for getPostCombsI
  mockr::local_mock(getPostCombsI = function(posterior_i) {
    return(0)
  })
  
  # Define a mock function for getPostProb
  mockr::local_mock(getPostProb = function(x, post_combs_i) {
    return(x)
  })
  
  # Define inputs
  posterior_i = 0
  contr_mat = list(contMat = matrix(c(0, 1), nrow = 1))
  crit_prob = 0.5
  
  # Call the function
  result = BayesMCPi(posterior_i, contr_mat, crit_prob)
  # Check the results
  expect_equal(result[["sign"]], 1)
  
  # Define inputs
  contr_mat = list(contMat = matrix(c(0, 0), nrow = 1))
  # Call the function
  result = BayesMCPi(posterior_i, contr_mat, crit_prob)
  # Check the results
  expect_equal(result[["sign"]], 0)

  # Define inputs
  contr_mat = list(contMat = matrix(c(0, 0.5), nrow = 1))
  # Call the function
  result = BayesMCPi(posterior_i, contr_mat, crit_prob)
  # Check the results
  expect_equal(result[["sign"]], 0)
})

#########################
# Tests for getPostProb #
#########################

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