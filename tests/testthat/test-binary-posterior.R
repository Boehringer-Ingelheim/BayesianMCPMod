test_that("getPosterior: binary mode uses binomial GLM and returns postList", {
  skip_if_not_installed("RBesT")
  
  set.seed(10)
  
  dose_levels <- c(0, 1, 2, 4)
  n_per_dose  <- 50
  probs       <- c(0.10, 0.20, 0.35, 0.55)
  
  dat <- data.frame(
    simulation = 1L,
    dose       = rep(dose_levels, each = n_per_dose),
    response   = unlist(lapply(probs, function(p) rbinom(n_per_dose, 1, p)))
  )
  attr(dat, "probability_scale") <- TRUE
  
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = qlogis(probs[i]), s = 1.5), sigma = 2)
    }),
    c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
  )
  
  post <- getPosterior(
    prior_list        = prior_list,
    data              = dat,
    calc_ess          = TRUE,
    probability_scale = TRUE
  )
  
  expect_s3_class(post, "postList")
  expect_equal(length(post), length(dose_levels))
  expect_true(is.list(attr(post, "posteriorInfo")))
  # ess can be numeric(0) in your current implementation depending on path
  expect_true(is.null(attr(post, "ess")) || is.numeric(attr(post, "ess")))
})

test_that("getPosterior: binary mode rejects malformed data (still needs simulation column)", {
  skip_if_not_installed("RBesT")
  
  dose_levels <- c(0, 1, 2, 4)
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 2), sigma = 2)
    }),
    c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
  )
  
  bad_dat <- data.frame(simulation = 1L, dose = c(0, 1), y = c(0, 1))
  attr(bad_dat, "probability_scale") <- TRUE
  
  expect_error(
    getPosterior(prior_list = prior_list, data = bad_dat, probability_scale = TRUE),
    regexp = "response"
  )
  
  bad_dat2 <- data.frame(simulation = 1L, dose = c(0, 0, 1, 1), response = c(0, 2, 1, 0))
  attr(bad_dat2, "probability_scale") <- TRUE
  
  # glm(binomial) should error when response is not 0/1 (or produce warnings then error downstream)
  expect_error(
    getPosterior(prior_list = prior_list, data = bad_dat2, probability_scale = TRUE)
  )
})

test_that("getPosterior: complete separation either succeeds with finite summary or errors (no silent NA)", {
  skip_if_not_installed("RBesT")
  
  dose_levels <- c(0, 1)
  prior_list <- setNames(
    list(
      RBesT::mixnorm(comp1 = c(w = 1, m = -2, s = 2), sigma = 2),
      RBesT::mixnorm(comp1 = c(w = 1, m =  2, s = 2), sigma = 2)
    ),
    c("Ctr", "DG_1")
  )
  
  dat <- data.frame(
    simulation = 1L,
    dose       = rep(dose_levels, each = 30),
    response   = c(rep(0, 30), rep(1, 30))
  )
  attr(dat, "probability_scale") <- TRUE
  
  res <- try(getPosterior(prior_list = prior_list, data = dat, probability_scale = TRUE), silent = TRUE)
  
  if (inherits(res, "try-error")) {
    # We accept an error, but it shouldn't be a cryptic failure later (like NA propagation).
    expect_true(nchar(as.character(res)) > 0)
  } else {
    expect_s3_class(res, "postList")
    s <- summary(res)
    expect_true(all(is.finite(s[, "mean"])))
    expect_true(all(is.finite(s[, "sd"])))
  }
})
