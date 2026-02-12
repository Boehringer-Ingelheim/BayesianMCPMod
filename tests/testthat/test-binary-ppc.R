test_that("Binary PPC: observed counts not extreme under posterior predictive (normal approx on logit)", {
  skip_if_not_installed("RBesT")
  
  set.seed(311)
  
  dose_levels <- c(0, 1, 2, 4)
  p_true      <- c(0.10, 0.20, 0.35, 0.55)
  n_per_dose  <- 120
  
  dat <- data.frame(
    simulation = 1L,
    dose       = rep(dose_levels, each = n_per_dose),
    response   = unlist(lapply(p_true, function(p) rbinom(n_per_dose, 1, p)))
  )
  attr(dat, "probability_scale") <- TRUE
  
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = qlogis(p_true[i]), s = 2.0), sigma = 2)
    }),
    c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
  )
  
  post <- getPosterior(prior_list = prior_list, data = dat, probability_scale = TRUE)
  
  y_obs <- as.numeric(tapply(dat$response, dat$dose, sum))
  
  smry <- summary(post)
  theta_mean <- as.numeric(smry[, "mean"])
  theta_sd   <- as.numeric(smry[, "sd"])
  
  expect_true(all(is.finite(theta_mean)))
  expect_true(all(is.finite(theta_sd)))
  expect_true(all(theta_sd > 0))
  
  ndraw <- 4000
  ppp <- numeric(length(dose_levels))
  
  for (j in seq_along(dose_levels)) {
    theta_draw <- rnorm(ndraw, mean = theta_mean[j], sd = theta_sd[j])
    p_draw <- plogis(theta_draw)
    y_rep <- rbinom(ndraw, size = n_per_dose, prob = p_draw)
    
    lo_tail <- mean(y_rep <= y_obs[j])
    hi_tail <- mean(y_rep >= y_obs[j])
    ppp[j] <- 2 * min(lo_tail, hi_tail)
  }
  
  # Not "essentially impossible under own posterior predictive"
  expect_true(all(ppp > 1e-4))
  # At least one dose looks reasonably typical
  expect_true(any(ppp > 0.02))
})
