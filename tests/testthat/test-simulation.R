# Import testthat unit testing framework
library(testthat) 

# Import functions to test
source("R/BMCPMod.R")

# Tests for assessDesign()

test_that("assessDesign() runs without errors", {
  # Define valid inputs
  n_patients <- c(10, 20, 30)
  mods <- getModel(c("linear", "quadratic")) 
  prior_list <- makePrior(mods, c(1,2,3))
  
  # Run function
  expect_error(assessDesign(n_patients, mods, prior_list), NA) 
})

test_that("assessDesign() output has expected structure", {
  # Define inputs
  n_patients <- c(10, 20, 30)
  mods <- getModel(c("linear", "quadratic"))
  prior_list <- makePrior(mods, c(1,2,3))  

  # Run function  
  result <- assessDesign(n_patients, mods, prior_list)

  # Validate output
  expect_is(result, "list")
  expect_named(result, c("linear", "quadratic"))
  expect_is(result$linear, "list")
  expect_is(result$quadratic, "list")
})

# Tests for getContrMat()

test_that("getContrMat() runs without errors", {
  # Define valid inputs
  mods <- getModel(c("linear", "quadratic"))
  dose_levels <- c(1,2,3)
  dose_weights <- c(10,20,30)
  prior_list <- makePrior(mods, dose_levels)

  # Run function
  expect_error(getContrMat(mods, dose_levels, dose_weights, prior_list), NA)
})

test_that("getContrMat() output has expected structure", {
  # Define inputs 
  mods <- getModel(c("linear", "quadratic"))
  dose_levels <- c(1,2,3)
  dose_weights <- c(10,20,30)
  prior_list <- makePrior(mods, dose_levels)

  # Run function
  result <- getContrMat(mods, dose_levels, dose_weights, prior_list)

  # Validate output
  expect_is(result, "list")
  expect_named(result, c("contMat", "muMat", "CorrMat"))
})

# Tests for getCritProb()

test_that("getCritProb() runs without errors", {
  # Define valid inputs
  mods <- getModel(c("linear", "quadratic"))
  dose_levels <- c(1,2,3)
  dose_weights <- c(10,20,30)

  # Run function
  expect_error(getCritProb(mods, dose_levels, dose_weights), NA)  
})

test_that("getCritProb() output is a numeric", {
  # Define inputs
  mods <- getModel(c("linear", "quadratic"))
  dose_levels <- c(1,2,3)
  dose_weights <- c(10,20,30)

  # Run function
  result <- getCritProb(mods, dose_levels, dose_weights)

  # Validate output
  expect_type(result, "double")
})
