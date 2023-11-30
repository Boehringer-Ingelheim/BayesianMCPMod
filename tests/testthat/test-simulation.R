# Test for getModelData
test_that("getModelData returns a data frame with correct columns", {

  mock_sim_data <- data.frame(
    simulation = 1:5,
    dose = rnorm(5),
    model1 = rnorm(5),
    model2 = rnorm(5)
  )
  mock_model_name <- "model1"

  result <- getModelData(sim_data = mock_sim_data, model_name = mock_model_name)
  
  # Assert that the result is a data frame
  expect_true(is.data.frame(result))
  # Assert that the data frame has the correct number of columns
  expect_equal(ncol(result), 3)
  # Assert that the data frame has the correct column names
  expect_equal(colnames(result), c("simulation", "dose", "response"))
})