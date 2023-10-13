test_that("getPosterior works correctly", {
  # Prepare test data and parameters
  data <- data.frame(simulation = rep(1, 4),
                     dose = c(0, 1, 2, 3),
                     response = c(10, 20, 30, 40))
  prior_list <- list(1, 2, 3, 4)
  mu_hat <- c(10, 20, 30, 40)
  sd_hat <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 1)
  
  # Test getPosterior function
  posterior_list <- getPosterior(data, prior_list, mu_hat, sd_hat)
  expect_type(posterior_list, "character")
  expect_s3_class(posterior_list, "postList")
})

test_that("getPosteriorI works correctly", {
  # Prepare test data and parameters
  data_i <- data.frame(dose = c(0, 1, 2, 3),
                       response = c(10, 20, 30, 40))
  prior_list <- list(1, 2, 3, 4)
  mu_hat <- c(10, 20, 30, 40)
  sd_hat <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 1)
  
  # Test getPosteriorI function
  post_list <- getPosteriorI(data_i, prior_list, mu_hat, sd_hat)
  expect_type(post_list, "character")
  expect_s3_class(post_list, "postList")
})

test_that("summary.postList works correctly", {
  # Prepare test data
  post_list <- list(
    Ctr = matrix(c(0.25, 10, 1), nrow = 3, ncol = 1),
    DG_1 = matrix(c(0.25, 20, 2), nrow = 3, ncol = 1),
    DG_2 = matrix(c(0.25, 30, 3), nrow = 3, ncol = 1),
    DG_3 = matrix(c(0.25, 40, 4), nrow = 3, ncol = 1)
  )
  class(post_list) <- "postList"
  
  # Test summary.postList function
  summary_tab <- summary.postList(post_list)
  expect_type(summary_tab, "character")
})
