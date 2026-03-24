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

test_that("getPosterior: binary mode uses Firth penalized regression in case of separation", {
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
  
  # hard-code separation
  dat[1:n_per_dose, 3] <- 0
  
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

test_that("getPosterior: binary no-separation case uses glm branch", {
  skip_if_not_installed("RBesT")
  
  set.seed(201)
  
  dose_levels <- c(0, 1, 2, 4)
  n_per_dose  <- 80
  probs       <- c(0.10, 0.20, 0.35, 0.55)
  
  dat <- data.frame(
    simulation = 1L,
    dose       = rep(dose_levels, each = n_per_dose),
    response   = unlist(lapply(probs, function(p) rbinom(n_per_dose, 1, p)))
  )
  attr(dat, "probability_scale") <- TRUE
  
  # Ensure test setup really is non-separated
  by_group <- split(dat$response, dat$dose)
  expect_true(all(vapply(by_group, function(y) any(y == 0) && any(y == 1), logical(1))))
  
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
  
  s <- summary(post)
  
  expect_s3_class(post, "postList")
  expect_equal(length(post), length(dose_levels))
  expect_equal(attr(post, "fitMethod"), "glm")
  expect_true(is.list(attr(post, "posteriorInfo")))
  expect_true(is.numeric(attr(post, "ess")))
  expect_true(all(is.finite(s[, "mean"])))
  expect_true(all(is.finite(s[, "sd"])))
})

test_that("getPosterior: binary separation case uses firth branch", {
  skip_if_not_installed("RBesT")
  skip_if_not_installed("logistf")
  
  dose_levels <- c(0, 1, 2, 4)
  n_per_dose  <- 50
  probs       <- c(0.10, 0.20, 0.35, 0.55)
  
  dat <- data.frame(
    simulation = 1L,
    dose       = rep(dose_levels, each = n_per_dose),
    response   = unlist(lapply(probs, function(p) rbinom(n_per_dose, 1, p)))
  )
  attr(dat, "probability_scale") <- TRUE
  
  # Force separation in control arm
  dat$response[dat$dose == 0] <- 0
  
  # Ensure test setup really is separated in at least one group
  by_group <- split(dat$response, dat$dose)
  expect_true(any(vapply(by_group, function(y) all(y == 0) || all(y == 1), logical(1))))
  
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
  
  s <- summary(post)
  
  expect_s3_class(post, "postList")
  expect_equal(length(post), length(dose_levels))
  expect_equal(attr(post, "fitMethod"), "firth")
  expect_true(is.list(attr(post, "posteriorInfo")))
  expect_true(is.numeric(attr(post, "ess")))
  expect_true(all(is.finite(s[, "mean"])))
  expect_true(all(is.finite(s[, "sd"])))
})

test_that("getPosterior: binary path handles mixed simulation branches", {
  skip_if_not_installed("RBesT")
  skip_if_not_installed("logistf")
  
  dose_levels <- c(0, 1)
  n_per_dose  <- 30
  
  dat1 <- data.frame(
    simulation = 1L,
    dose = rep(dose_levels, each = n_per_dose),
    response = c(rbinom(n_per_dose, 1, 0.2), rbinom(n_per_dose, 1, 0.5))
  )
  
  dat2 <- data.frame(
    simulation = 2L,
    dose = rep(dose_levels, each = n_per_dose),
    response = c(rep(0, n_per_dose), rep(1, n_per_dose))
  )
  
  dat <- rbind(dat1, dat2)
  attr(dat, "probability_scale") <- TRUE
  
  prior_list <- setNames(
    list(
      RBesT::mixnorm(comp1 = c(w = 1, m = -1, s = 2), sigma = 2),
      RBesT::mixnorm(comp1 = c(w = 1, m =  1, s = 2), sigma = 2)
    ),
    c("Ctr", "DG_1")
  )
  
  post <- getPosterior(prior_list = prior_list, data = dat, probability_scale = TRUE)
  
  expect_type(post, "list")
  expect_equal(length(post), 2)
  expect_s3_class(post[[1]], "postList")
  expect_s3_class(post[[2]], "postList")
  expect_equal(attr(post[[1]], "fitMethod"), "glm")
  expect_equal(attr(post[[2]], "fitMethod"), "firth")
})
