library(testthat) 
library(mockery)

test_that("simulateData returns expected output", {

  # Set up
  n_patients <- c(5, 10, 15)
  dose_levels <- c(1, 5, 10) 
  sd <- 2
  mods <- getTestMods() # Function to generate test mods object

  # Exercise  
  result <- simulateData(n_patients, dose_levels, sd, mods)

  # Verify
  expect_s3_class(result, "data.frame")
  expect_named(result, c("simulation", "ptno", "dose", "response")) 
  expect_equal(nrow(result), sum(n_patients))
  expect_true(all(result$dose %in% dose_levels))

})

test_that("getModelData returns expected model data", {

  # Set up
  sim_data <- getSimData() # Function to generate test data
  model_name <- "linear"

  # Exercise
  result <- getModelData(sim_data, model_name)

  # Verify
  expect_s3_class(result, "data.frame")
  expect_named(result, c("simulation", "dose", "response"))
  expect_equal(ncol(result), 3)
  expect_true(all(result$response == sim_data[[model_name]]))

})
