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

# Test predictModelFit function
test_that("predictModelFit works correctly", {
  model_fit <- list(
    model = "emax",
    coeffs = c(e0 = 0, eMax = 1, ed50 = 2),
    dose_levels = c(0, 1, 2, 3, 4, 5)
  )
  
  pred_vals <- predictModelFit(model_fit)
  expect_type(pred_vals, "double")
})
