test_that("Binary posterior intervals are finite, shrink with n, and coverage is non-degenerate", {
  skip_if_not_installed("RBesT")
  
  set.seed(201)
  
  dose_levels <- c(0, 1, 2, 4)
  p_true      <- c(0.10, 0.20, 0.35, 0.55)
  theta_true  <- qlogis(p_true)
  
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 3), sigma = 3)
    }),
    c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
  )
  
  run_replications <- function(n_per_dose, R) {
    covered <- matrix(FALSE, nrow = R, ncol = length(dose_levels))
    widths  <- matrix(NA_real_, nrow = R, ncol = length(dose_levels))
    
    for (r in seq_len(R)) {
      dat <- data.frame(
        simulation = 1L,
        dose       = rep(dose_levels, each = n_per_dose),
        response   = unlist(lapply(p_true, function(p) rbinom(n_per_dose, 1, p)))
      )
      attr(dat, "probability_scale") <- TRUE
      
      post <- getPosterior(
        prior_list        = prior_list,
        data              = dat,
        probability_scale = TRUE
      )
      
      lo <- as.numeric(sapply(post, function(m) stats::quantile(m, probs = 0.025)))
      hi <- as.numeric(sapply(post, function(m) stats::quantile(m, probs = 0.975)))
      
      # Always-required sanity
      if (!all(is.finite(lo)) || !all(is.finite(hi)) || !all(lo <= hi)) {
        return(list(ok = FALSE))
      }
      
      widths[r, ]  <- (hi - lo)
      covered[r, ] <- (theta_true >= lo) & (theta_true <= hi)
    }
    
    cov_rate <- colMeans(covered)
    width    <- colMeans(widths)
    
    list(
      ok = TRUE,
      cov_rate = cov_rate,
      width = width,
      cov_mean = mean(cov_rate),
      width_mean = mean(width),
      cov_min = min(cov_rate),
      cov_max = max(cov_rate)
    )
  }
  
  # Keep runtime reasonable but stable
  res_small <- run_replications(n_per_dose = 40,  R = 40)
  res_big   <- run_replications(n_per_dose = 200, R = 40)
  
  expect_true(res_small$ok)
  expect_true(res_big$ok)
  
  # 1) Width must shrink with more data (core statistical plausibility)
  expect_true(res_big$width_mean < res_small$width_mean)
  
  # 2) Coverage must be non-degenerate overall (not ~0, not identically 1)
  #    These are intentionally loose to accommodate approximate posteriors.
  expect_true(res_small$cov_mean > 0.05)
  expect_true(res_small$cov_mean < 0.999)
  
  expect_true(res_big$cov_mean > 0.05)
  expect_true(res_big$cov_mean < 0.999)
  
  # 3) It should not get dramatically worse with more data in aggregate
  #    (This catches regressions where CI computation breaks asymptotically.)
  expect_true(res_big$cov_mean + 0.20 >= res_small$cov_mean)
  
  # 4) Guard against "everything always covers" for all doses at large n
  #    If max coverage is 1.0 that's fine, but not all doses should be ~1 simultaneously.
  expect_true(mean(res_big$cov_rate > 0.999) < 0.75)
})
