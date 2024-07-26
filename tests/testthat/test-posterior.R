# Function to test if covariance matrix is symmerical
test_covmatrix_symmetry <- function (

    posterior,
    mu_hat

  ) {

  lapply(1:(length(mu_hat)-1), function(x) {

    expect_equal(posterior$covmat$Comp1[row(posterior$covmat$Comp1) + x == col(posterior$covmat$Comp1)],
                 posterior$covmat$Comp1[row(posterior$covmat$Comp1)     == col(posterior$covmat$Comp1) + x])

    expect_equal(length(posterior$covmat$Comp1[row(posterior$covmat$Comp1) + x == col(posterior$covmat$Comp1)]),
                 length(posterior$covmat$Comp1[row(posterior$covmat$Comp1) + x       == col(posterior$covmat$Comp1)]))

  })

}

test_that("getPosterior works correctly", {
  dummy_data <- getModelData(data, names(mods)[1])
  prior_list_noRBesT <- list(1,2,3,4)

  expect_error(getPosterior(
    prior_list = prior_list_noRBesT,
    mu_hat = mu_hat,
    S_hat = se_hat_matrix,
    calc_ess = FALSE
  ))

  # Test getPosterior function
  posterior_list <- getPosterior(
    prior_list = prior_list_matrix,
    mu_hat = mu_hat,
    S_hat = se_hat_vector,
    calc_ess = FALSE
  )
  expect_type(posterior_list, "list")
  expect_s3_class(posterior_list, "postList")
})

test_that("getPriorList input parameters do work as intented", {
  set.seed(8080)
  dataset     <- dplyr::filter(testdata, bname == "BRINTELLIX")
  histcontrol <- dplyr::filter(dataset, dose == 0, primtime == 8, indication == "MAJOR DEPRESSIVE DISORDER",protid!=6)

  ##Create MAP Prior
  hist_data <- data.frame(
    trial = histcontrol$nctno,
    est   = histcontrol$rslt,
    se    = histcontrol$se,
    sd    = histcontrol$sd,
    n     = histcontrol$sampsize)

  dose_levels <- c(0, 2.5, 5, 10)

  # call without parameter
  expect_error(getPriorList())
  # only passing the historical data
  expect_error(getPriorList(hist_data = hist_data))
  # passing both needed parameters
  expect_type(getPriorList(
    hist_data = hist_data,
    dose_levels = dose_levels,
    robust_weight = 0.5
  ), "list")
  # passing wrong format for hist_data
  expect_error(getPriorList(
    hist_data = testdata,
    dose_levels = dose_levels,
    robust_weight = 0.5
  ))
  # passing wrong format for dose_levels
  expect_error(getPriorList(
    hist_data = hist_data,
    dose_levels = c("hello", "world"),
    robust_weight = 0.5
  ))

})

test_that("getPosteriorI works correctly", {
  # Prepare test data and parameters
  data_i <- data.frame(
    dose = c(0, 1, 2, 3, 4),
    response = c(10, 20, 30, 40, 50)
  )

  #prior_list <- list(1, 2, 3, 4)
  #mu_hat <- c(10, 20, 30, 40)
  #se_hat <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 1)

  # Test getPosteriorI function
  post_list <- getPosteriorI(data_i = data_i, prior_list = prior_list_matrix, mu_hat = mu_hat, se_hat = se_hat_vector)
  expect_type(post_list, "list")
  expect_s3_class(post_list, "postList")

  # Test mu_hat and sd_hat both null branch
  post_list <- getPosteriorI(data_i, prior_list_matrix, NULL, NULL)
  expect_type(post_list, "list")
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

test_that("getPostCombsI returns an object with correct attributes", {
  posterior_i <- list(
    matrix(c(2, 1, 2), nrow = 3),
    matrix(c(2, 1, 2), nrow = 3)
  )
  result <- getPostCombsI(posterior_i)

  expect_true(is.list(result))
  expect_equal(length(result), 3)
  expect_equal(names(result), c("weights", "means", "vars"))
  expect_equal(result$weights, 4)
  expect_true(all(result$vars == c(4, 4)))
})

test_that("priorList2priorMix works correctly", {

  # function call without parameters
  expect_error(priorList2priorMix())

  # test priorList2priorMix function
  prior_mix <- priorList2priorMix(prior_list_matrix)
  expect_type(prior_mix, "list")
  expect_length(prior_mix, 3)

})

test_that("postMix2posteriorList works correctly", {

  # create posterior with matrix
  prior_mix <- priorList2priorMix(prior_list_matrix)

  posterior <- DoseFinding::mvpostmix(
    priormix = prior_mix,
    mu_hat   = mu_hat,
    S_hat    = se_hat_matrix)

  # create posterior with vector
  posterior_vector <- getPosterior(
    prior_list = prior_list_matrix,
    mu_hat = mu_hat,
    S_hat = se_hat_vector,
    calc_ess = FALSE
  )

  ### test postMix2posteriorList function - se_hat_matrix only zeros
  posterior_list <- postMix2posteriorList(posterior, prior_list_matrix, calc_ess = FALSE)
  expect_type(posterior_list, "list")
  expect_s3_class(posterior_list, "postList")

  test_covmatrix_symmetry(posterior, mu_hat)

  # compare posterior result object with matrix to object with vector
  expect_length(posterior_list, length(posterior_vector))
  expect_length(attr(posterior_list, "covariance matrices"), length(attr(posterior_vector, "full covariance matrices")))
  expect_type(posterior_list, typeof(posterior_vector))
  expect_s3_class(posterior_list, S3Class(posterior_vector))


  ### test postMix2posteriorList function - se_hat_matrix not only zeros
  mvpostmix_noZero <- DoseFinding::mvpostmix(
    priormix = prior_mix,
    mu_hat   = mu_hat,
    S_hat    = se_hat_matrix2)

  posterior_noZero <- getPosterior(
  prior_list = prior_list_matrix,
  mu_hat     = mu_hat,
  S_hat      = se_hat_matrix2,
  calc_ess   = FALSE
  )

  expect_type(posterior_noZero, "list")
  expect_s3_class(posterior_noZero, "postList")

  test_covmatrix_symmetry(mvpostmix_noZero, mu_hat)

  # compare posterior result object with matrix to object with vector
  expect_length(posterior_noZero, length(posterior_vector))
  expect_length(attr(posterior_noZero, "covariance matrices"), length(attr(posterior_vector, "full covariance matrices")))
  expect_type(posterior_noZero, typeof(posterior_vector))
  expect_s3_class(posterior_noZero, S3Class(posterior_vector))


  ### test similarity of results for posterior from a matrix compared to posterior from a vector with square roots
  se_hat <- c(1, 2, 3, 4, 5)
  S_hat  <- diag(se_hat)

  posterior_matrix_S <- getPosterior(
    prior_list = prior_list_matrix,
    mu_hat     = mu_hat,
    S_hat      = S_hat,
    calc_ess   = FALSE
  )

  posterior_vector_se <- getPosterior(
    prior_list = prior_list_matrix,
    mu_hat     = mu_hat,
    S_hat      = sqrt(se_hat),
    calc_ess   = FALSE
  )

  expect_equal(posterior_matrix_S$Ctr, posterior_vector_se$Ctr)
  expect_equal(posterior_matrix_S$DG_1, posterior_vector_se$DG_1)
  expect_equal(posterior_matrix_S$DG_2, posterior_vector_se$DG_2)
  expect_equal(posterior_matrix_S$DG_3, posterior_vector_se$DG_3)
  expect_equal(posterior_matrix_S$DG_4, posterior_vector_se$DG_4)




})
