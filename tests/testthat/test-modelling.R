# Test predictModelFit function
test_that("predictModelFit works correctly", {

  # Test data
  test_data <- data.frame(
    simulation = rep(1, 6),
    dose = c(0, 1, 2, 3, 4, 5),
    response = c(0, 1, 2, 3, 4, 5)
  )

  # Mock getPosterior function
  getPosterior <- function(data, prior_list, mu_hat, sd_hat) {
    list(
      means = c(0, 1, 2, 3, 4, 5),
      vars = c(1, 1, 1, 1, 1, 1),
      weights = c(1, 1, 1, 1, 1, 1)
    )
  }
  
  model_fit <- list(
    model = "emax",
    coeffs = c(e0 = 0, eMax = 1, ed50 = 2),
    dose_levels = c(0, 1, 2, 3, 4, 5)
  )
  
  pred_vals <- predictModelFit(model_fit)
  expect_type(pred_vals, "double")
})

# Test for addModelWeights
test_that("addModelWeights works correctly in a simple case", {
  # Create a test case
  test_model_fits1 <- list(
    model1 = list(gAIC = 1),
    model2 = list(gAIC = 1),
    model3 = list(gAIC = 1)
  )
  expected_output1 <- list(
    list(gAIC = 1, model_weight = 1/3),
    list(gAIC = 1, model_weight = 1/3),
    list(gAIC = 1, model_weight = 1/3)
  )
  
  # Call the function with the test case
  result <- addModelWeights(model_fits = test_model_fits1)
  # Assert the expected behavior
  expect_true(is.list(result))
  expect_equal(length(result), length(test_model_fits))
  expect_true(all(sapply(result, function(x) "model_weight" %in% names(x))))
  expect_equal(result, expected_output1)

  test_model_fits2 <- list(
    model1 = list(gAIC = 100),
    model2 = list(gAIC = 0),
    model3 = list(gAIC = 100)
  )
  expected_output2 <- list(
    list(gAIC = 100, model_weight = 0),
    list(gAIC = 0, model_weight = 1),
    list(gAIC = 100, model_weight = 0)
  )
  result <- addModelWeights(model_fits = test_model_fits2)
  expect_equal(result, expected_output2, tolerance = 1e-10)
})
