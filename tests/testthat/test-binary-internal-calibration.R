test_that("Binary internal calibration: posterior mean ranks with empirical logits; SD tracks binomial information", {
  skip_if_not_installed("RBesT")
  
  set.seed(310)
  
  dose_levels <- c(0, 1, 2, 4)
  
  # Use separated but not extreme probabilities to avoid separation artifacts
  p_true <- c(0.12, 0.20, 0.33, 0.52)
  theta_true <- qlogis(p_true)
  
  # Mild priors centered near truth to stabilize CI; still allows data to dominate
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = theta_true[i], s = 2.0), sigma = 2)
    }),
    c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
  )
  
  run_once <- function(n_per_dose) {
    dat <- data.frame(
      simulation = 1L,
      dose       = rep(dose_levels, each = n_per_dose),
      response   = unlist(lapply(p_true, function(p) rbinom(n_per_dose, 1, p)))
    )
    attr(dat, "probability_scale") <- TRUE
    
    post <- getPosterior(prior_list = prior_list, data = dat, probability_scale = TRUE)
    smry <- summary(post)
    
    # Empirical proportions/logits with continuity correction to avoid +/-Inf
    p_hat <- tapply(dat$response, dat$dose, mean)
    eps <- 1/(n_per_dose + 1)
    p_hat_cc <- pmin(pmax(as.numeric(p_hat), eps), 1 - eps)
    logit_hat <- qlogis(p_hat_cc)
    
    list(
      mean = as.numeric(smry[, "mean"]),
      sd   = as.numeric(smry[, "sd"]),
      logit_hat = as.numeric(logit_hat),
      p_hat = as.numeric(p_hat_cc)
    )
  }
  
  # Compare two sample sizes to test shrinkage + scaling
  small_n <- 40
  big_n   <- 200
  
  out_small <- run_once(small_n)
  out_big   <- run_once(big_n)
  
  # --- Means: plausible & monotone (ranking) ---
  # rank(posterior mean) should match rank(empirical logit)
  expect_equal(order(out_small$mean), order(out_small$logit_hat))
  expect_equal(order(out_big$mean),   order(out_big$logit_hat))
  
  # Means should be reasonably close to empirical logits (looser at small n)
  expect_equal(out_small$mean, out_small$logit_hat, tolerance = 0.60)
  expect_equal(out_big$mean,   out_big$logit_hat,   tolerance = 0.35)
  
  # --- SD sanity ---
  expect_true(all(is.finite(out_small$sd)))
  expect_true(all(is.finite(out_big$sd)))
  expect_true(all(out_small$sd > 0))
  expect_true(all(out_big$sd > 0))
  
  # --- Shrinkage with n (must hold in aggregate) ---
  expect_true(mean(out_big$sd) < mean(out_small$sd))
  
  # --- Binomial information scaling on logit scale ---
  # For Bernoulli with logit link, Fisher info ~ n p (1-p), so sd(theta) ~ 1/sqrt(n p(1-p))
  approx_sd <- function(n, p) 1 / sqrt(n * p * (1 - p))
  
  # Compare ratio of SDs between small and big n: should be ~sqrt(big/small) up to slack
  ratio_obs <- out_small$sd / out_big$sd
  ratio_exp <- sqrt(big_n / small_n) * approx_sd(small_n, out_small$p_hat) / approx_sd(big_n, out_big$p_hat)
  # ratio_exp simplifies close to sqrt(big/small) if p's match; keep robust by computing anyway.
  
  # Require positive association between observed and expected scaling across doses
  expect_true(stats::cor(ratio_obs, ratio_exp) > 0.3)
  
  # Also require median ratio in a plausible band around sqrt(big/small) (~2.236)
  target <- sqrt(big_n / small_n)
  expect_true(stats::median(ratio_obs) > 0.9 * target)
  expect_true(stats::median(ratio_obs) < 1.6 * target)
})
