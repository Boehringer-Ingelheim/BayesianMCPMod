##########################
# Tests for assessDesign #
##########################

test_that("base case input throws no error", {
  
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
  
  expect_equal(
    names(eval_design), names(mods)
  )
  
  expect_equal(
    attr(eval_design$linear$BayesianMCP, "dim")[1],
    n_sim
  )
  
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
