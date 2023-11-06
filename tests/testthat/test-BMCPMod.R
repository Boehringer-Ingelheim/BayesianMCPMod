# Load the necessary library
library(testthat)

# Create minimal test case
n_hist_trials = 1

hist_data <- data.frame(
  trial = seq(1, n_hist_trials, 1),
  est   = rep(-1, n_hist_trials),
  se    = rep(1, n_hist_trials),
  sd    = rep(1, n_hist_trials),
  n     = rep(1, n_hist_trials)
)

n_patients <- c(2, 2)
dose_levels <- c(0, 2.5)

mods <- DoseFinding::Mods(
  linear = NULL, 
  doses = dose_levels
)

prior_list <- getPriorList(
  hist_data   = hist_data,
  dose_levels = dose_levels,
  robustify_weight = 0.5
)

# Tests for n_patients parameter
test_that("assessDesign correctly validates n_patients parameter input", {

  expect_error(
    assessDesign(
      n_patients = NULL, 
      mods = mods, 
      prior_list = prior_list
    ),
    "not be NULL",
    ignore.case = T
  )
  
  expect_error(
    assessDesign(
      n_patients = n_patients[-1], 
      mods = mods, 
      prior_list = prior_list
    ),
    "match number of dose levels",
    ignore.case = T
  )
  
  expect_error(
    assessDesign(
      n_patients = rep(1, length(n_patients)),
      mods = mods,
      prior_list = prior_list
    ),
    "at least one value must be greater than 1",
    ignore.case = T
  )
  
  expect_no_error(
    assessDesign(
      n_patients = n_patients - c(0, 1),
      mods = mods,
      prior_list = prior_list
    )
  )

  expect_no_error(
    assessDesign(
      n_patients = n_patients, 
      mods = mods, 
      prior_list = prior_list
    )
  )
})
# 
# # Test for mods parameter
# test_that("mods parameter works correctly", {
#   # Call the function with the test case
#   result <- assessDesign(
#     n_patients = n_patients, 
#     mods = mods, 
#     prior_list = prior_list
#   )
#   # Assert the expected behavior
#   expect_equal(class(result$mods), "Mods")
# })
# 
# # Test for prior_list parameter
# test_that("prior_list parameter works correctly", {
#   # Create a test case
#   test_prior_list <- list()
#   # Call the function with the test case
#   result <- assessDesign(n_patients = NULL, mods = NULL, prior_list = test_prior_list)
#   # Assert the expected behavior
#   expect_equal(class(result$prior_list), "list")
# })