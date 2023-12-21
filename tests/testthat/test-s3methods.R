test_that("print.BayesianMCPMod works as intented", {
  expect_error(print.BayesianMCPMod())
  
})


test_that("print.BayesianMCP works as intented", {
  expect_error(print.BayesianMCP())
  
  
  
})

test_that("predict.ModelFits works as intented", {
  expect_error(predict.ModelFits())
  
})

test_that("s3 postList functions work as intented", {
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
  post_test_list <- getPriorList(
    hist_data = hist_data,
    dose_levels = dose_levels,
    robust_weight = 0.5)
  expect_error(summary.postList())
  expect_type(summary.postList(post_test_list), "double")
  expect_error(print.postList())
  expect_type(print(post_test_list), "list")
  expect_type(print.postList(post_test_list), "list")
  # expect_true(names(print(post_test_list)) == c("Summary of Posterior Distributions",
  #                                               "Maximum Difference to Control and Dose Group",
  #                                               "Posterior Distributions"))
})

test_that("test modelFits s3 methods", {
  
  model_shapes <- colnames(contr_mat$contMat)
  dose_levels  <- as.numeric(rownames(contr_mat$contMat))
  
  model_fits  <- getModelFits(
    models      = model_shapes,
    dose_levels = dose_levels,
    posterior   = posterior_list,
    simple      = simple)
  
  pred <- predict(model_fits)
  pred_dosage <- predict(model_fits, doses = dose_levels)
  
  expect_type(pred, "list")
  expect_true(is.null(attr(pred, "doses")))
  expect_identical(attr(pred_dosage, "doses"), dose_levels)
  expect_type(print(model_fits), "double")
})
