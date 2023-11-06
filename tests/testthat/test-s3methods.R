# Import testthat
library(testthat) 

# Import functions 
source("R/s3methods.R")

test_that("print.BayesianMCPMod prints model summary", {
  # Arrange
  x <- create_BayesianMCPMod_object()
  
  # Act
  print_output <- capture_print(print.BayesianMCPMod(x))

  # Assert 
  expect_match(print_output, "Bayesian Multiple Comparison Procedure", fixed=TRUE)
  expect_match(print_output, "Model Significance Frequencies", fixed=TRUE)
})

test_that("print.BayesianMCP prints power and n_sim", {
  # Arrange
  x <- create_BayesianMCP_object()
  
  # Act
  print_output <- capture_print(print.BayesianMCP(x))

  # Assert
  expect_match(print_output, "Estimated Success Rate:", fixed=TRUE)
  expect_match(print_output, "N Simulations:", fixed=TRUE)  
})
