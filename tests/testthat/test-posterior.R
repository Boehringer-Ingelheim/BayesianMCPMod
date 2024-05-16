test_that("getPosterior works correctly", {
  dummy_data <- getModelData(data, names(mods)[1])
  
  # Test getPosterior function
  posterior_list <- getPosterior(
    data = getModelData(data, names(mods)[1]),
    prior_list = prior_list
  )
  expect_type(posterior_list, "list")
  expect_s3_class(posterior_list, "postList")
})

test_that("getPriorList input parameters do work as intented", {
  # setup
  library(clinDR)
  library(dplyr)
  set.seed(8080)
  data("metaData")
  testdata    <- as.data.frame(metaData)
  dataset     <- filter(testdata, bname == "BRINTELLIX")
  histcontrol <- filter(dataset, dose == 0, primtime == 8, indication == "MAJOR DEPRESSIVE DISORDER",protid!=6)
  
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
    dose = c(0, 1, 2, 3),
    response = c(10, 20, 30, 40)
  )
  
  prior_list <- list(1, 2, 3, 4)
  mu_hat <- c(10, 20, 30, 40)
  se_hat <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 1)
  
  # Test getPosteriorI function
  post_list <- getPosteriorI(data_i, prior_list, mu_hat, se_hat)
  expect_type(post_list, "list")
  expect_s3_class(post_list, "postList")
  
  # Test mu_hat and sd_hat both null branch
  post_list <- getPosteriorI(data_i, prior_list, NULL, NULL)
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

test_that("getPriorMix works correctly", {
  
  # function call without parameters
  expect_error(createPriorMix())
  
  # test createPriorMix function
  prior_mix <- createPriorMix(prior_list_covmat)
  expect_type(prior_mix, "list")
  expect_length(prior_mix, 3)
  
})

test_that("mvpostmix works correctly", {
  
  prior_mix <- createPriorMix(prior_list_covmat)
  
  # test mvpostmix function
  expect_error(mvpostmix(prior_mix, mu_hat, se_hat_vector))
  expect_no_error(mvpostmix(prior_mix, mu_hat, se_hat_matrix))
  posterior <- mvpostmix(prior_mix, mu_hat, se_hat_matrix)
  expect_type(posterior, "list")
  expect_length(posterior, 3)
  
})

test_that("getPosteriorOutput works correctly", {
  
  # create posterior with matrix
  prior_mix <- createPriorMix(prior_list_covmat)
  
  posterior <- DoseFinding::mvpostmix(
    priormix = prior_mix,
    mu_hat   = mu_hat,
    S_hat    = se_hat_matrix)
  
  # create posterior with vector
  posterior_vector <- getPosteriorI(
    prior_list = prior_list_covmat,
    mu_hat = mu_hat,
    se_hat = se_hat_vector,
    calc_ess = FALSE
  )
  
  # test getPosteriorOutput function
  posterior_list <- postmix2RBesT(posterior, prior_list_covmat, calc_ess = FALSE)
  expect_type(posterior_list, "list")
  expect_s3_class(posterior_list, "postList")
  
  lapply(1:(length(dose_levels_covmat)-1), function(x) {
    
    expect_equal(posterior$covmat$Comp1[row(posterior$covmat$Comp1) + x == col(posterior$covmat$Comp1)], 
                 posterior$covmat$Comp1[row(posterior$covmat$Comp1)     == col(posterior$covmat$Comp1) + x])
    
    expect_equal(length(which(posterior$covmat$Comp1[row(posterior$covmat$Comp1) + x == col(posterior$covmat$Comp1)] >= 0)),
                 length(posterior$covmat$Comp1[row(posterior$covmat$Comp1) + x       == col(posterior$covmat$Comp1)]))
    
  })
  
  lapply(seq_along(prior_list), function(i) {
    
    sapply(seq_along(posterior$weights), function(j) {
        
      expect_in(prior_list[[i]]["m",], prior_mix[[2]][[j]][i])
        
    })
      
  })
  
  expect_equal(length(posterior$weights), prod(lengths(prior_list_covmat)/nrow(prior_list_covmat$Ctr)))
  
  # compare posterior result object with matrix to object with vector
  expect_length(posterior_list, length(posterior_vector))
  expect_type(posterior_list, typeof(posterior_vector))
  expect_s3_class(posterior_list, S3Class(posterior_vector))
  
})