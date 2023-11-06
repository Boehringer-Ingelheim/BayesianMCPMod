# Import functions to test 
source("R/BMCPMod.R")

# Import testthat for unit testing
library(testthat) 

# Tests for assessDesign
test_that("assessDesign returns correct class", {
  expect_s3_class(assessDesign(n_patients, mods, prior_list), "BayesianMCPMod") 
})

test_that("assessDesign handles errors", {
  expect_error(assessDesign(n_patients, mods, prior_list, n_sim = "string"))
  expect_error(assessDesign(n_patients, mods, prior_list, alpha_crit_val = -1))
})

# Tests for getContrMat
test_that("getContrMat returns correct class", {
  expect_s3_class(getContrMat(mods, dose_levels, dose_weights, prior_list), "optContr")
})

test_that("getContrMat handles errors", {
  expect_error(getContrMat(mods, dose_levels, dose_weights, prior_list = 1))
})

# Tests for getCritProb
test_that("getCritProb returns a probability", {
  expect_gte(getCritProb(mods, dose_levels, dose_weights), 0)
  expect_lte(getCritProb(mods, dose_levels, dose_weights), 1)  
})

test_that("getCritProb handles errors", {
  expect_error(getCritProb(mods, dose_levels, dose_weights, alpha_crit_val = -1))
})

# Tests for performBayesianMCPMod
test_that("performBayesianMCPMod returns correct class", {
  expect_s3_class(performBayesianMCPMod(posteriors_list, contr_mat, crit_prob), "BayesianMCPMod")
})

test_that("performBayesianMCPMod handles errors", {
  expect_error(performBayesianMCPMod(posteriors_list, contr_mat, crit_prob = -1))
})

# Tests for addSignificance
test_that("addSignificance adds significant attribute", {
  model_fits <- list(a = 1, b = 2)
  expect_named(addSignificance(model_fits, c(TRUE, FALSE)), "significant") 
})

# Tests for performBayesianMCP
test_that("performBayesianMCP returns correct class", {
  expect_s3_class(performBayesianMCP(posteriors_list, contr_mat, crit_prob), "BayesianMCP")  
})

test_that("performBayesianMCP handles errors", {
  expect_error(performBayesianMCP(posteriors_list, contr_mat, crit_prob = -1))
})
