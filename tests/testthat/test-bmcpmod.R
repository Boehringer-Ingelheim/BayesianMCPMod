##########################
# Tests for assessDesign #
##########################

test_that("base case input throws no error", {
  expect_no_error(
    assessDesign(
      n_patients = n_patients, 
      mods = mods, 
      prior_list = prior_list
    )
  )
})


### n_patients parameter ###
test_that("assessDesign validates n_patients parameter input and give appropriate error messages", {
  
  expect_error(
    assessDesign(n_patients = NULL, mods = mods, prior_list = prior_list),
    "n_patients should not be NULL", ignore.case = T
  )
  
  expect_error(
    assessDesign(n_patients = list(), mods = mods, prior_list = prior_list),
    "n_patients should be a vector", ignore.case = T
  )
  
  expect_error(
    assessDesign(n_patients = c("2", "2"), mods = mods, prior_list = prior_list),
    "n_patients should be numeric", ignore.case = T
  )
  
  expect_error(
    assessDesign(n_patients = n_patients[-1], mods = mods, prior_list = prior_list),
    "length of n_patients should equal number of dose groups", ignore.case = T
  )
  
  expect_error(
    assessDesign(n_patients = rep(1, length(n_patients)), mods = mods, prior_list = prior_list),
    "at least one element in n_patients needs to be > 1", ignore.case = T
  )
})

### mods parameter ###
test_that("assessDesign validates mods parameter input and give appropriate error messages", {
  
  # assertions that aren't tested here for sake of brevity
  # mods should not be NULL
  # mods should be of class "Mods" from {DoseFinding}
  # length(n_patients) == length(attributes(mods)$doses) is commutative, so testing here is redundant
 
})