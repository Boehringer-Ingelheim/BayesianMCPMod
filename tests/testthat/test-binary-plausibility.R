test_that("Binary posterior: posterior means are close to empirical logits (vague prior, large n)", {
  skip_if_not_installed("RBesT")
  
  set.seed(101)
  
  dose_levels <- c(0, 1, 2, 4)
  n_per_dose  <- 400
  p_true      <- c(0.10, 0.20, 0.35, 0.55)
  
  dat <- data.frame(
    simulation = 1L,
    dose       = rep(dose_levels, each = n_per_dose),
    response   = unlist(lapply(p_true, function(p) rbinom(n_per_dose, 1, p)))
  )
  attr(dat, "probability_scale") <- TRUE
  
  # Vague priors on logit scale (so data dominates)
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 10), sigma = 10)
    }),
    c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
  )
  
  post <- getPosterior(prior_list = prior_list, data = dat, probability_scale = TRUE)
  smry <- summary(post)
  
  # empirical logits (avoid +/-Inf)
  p_hat <- tapply(dat$response, dat$dose, mean)
  eps <- 1/(n_per_dose + 1)  # gentle continuity correction
  p_hat_cc <- pmin(pmax(p_hat, eps), 1 - eps)
  logit_hat <- qlogis(p_hat_cc)
  
  # posterior means should be near empirical logits when n is large and prior vague
  # Tolerance: 0.25 on logit scale is fairly tight for Monte Carlo variability
  expect_equal(
    as.numeric(smry[, "mean"]),
    as.numeric(logit_hat),
    tolerance = 0.25
  )
  
  # sanity: posterior SDs should be finite and not huge with n=400
  expect_true(all(is.finite(smry[, "sd"])))
  expect_true(all(smry[, "sd"] < 2))
})


test_that("Binary posterior: uncertainty shrinks with more data (same p_true)", {
  skip_if_not_installed("RBesT")
  
  set.seed(102)
  
  dose_levels <- c(0, 1, 2, 4)
  p_true      <- c(0.12, 0.18, 0.30, 0.50)
  
  make_dat <- function(n_per_dose) {
    d <- data.frame(
      simulation = 1L,
      dose       = rep(dose_levels, each = n_per_dose),
      response   = unlist(lapply(p_true, function(p) rbinom(n_per_dose, 1, p)))
    )
    attr(d, "probability_scale") <- TRUE
    d
  }
  
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 10), sigma = 10)
    }),
    c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
  )
  
  post_small <- getPosterior(prior_list = prior_list, data = make_dat(40), probability_scale = TRUE)
  post_big   <- getPosterior(prior_list = prior_list, data = make_dat(400), probability_scale = TRUE)
  
  sd_small <- as.numeric(summary(post_small)[, "sd"])
  sd_big   <- as.numeric(summary(post_big)[, "sd"])
  
  # Big sample should have smaller SD in aggregate (not necessarily every single dose)
  expect_true(mean(sd_big) < mean(sd_small))
})


test_that("Binary modeling: fitted probabilities are plausible vs. generating truth (strong signal)", {
  skip_if_not_installed("RBesT")
  skip_if_not_installed("DoseFinding")
  
  set.seed(103)
  
  dose_levels <- c(0, 1, 2, 4)
  n_per_dose  <- 250
  p_true      <- c(0.08, 0.18, 0.35, 0.65)
  
  dat <- data.frame(
    simulation = 1L,
    dose       = rep(dose_levels, each = n_per_dose),
    response   = unlist(lapply(p_true, function(p) rbinom(n_per_dose, 1, p)))
  )
  attr(dat, "probability_scale") <- TRUE
  
  # Mildly informative priors near truth (keeps test stable)
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = qlogis(p_true[i]), s = 1.5), sigma = 2)
    }),
    c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
  )
  
  post <- getPosterior(prior_list = prior_list, data = dat, probability_scale = TRUE)
  
  mods <- DoseFinding::Mods(
    linear = NULL,
    emax   = c(0.5, 1.2),
    doses  = dose_levels,
    maxEff = 2
  )
  
  fits <- getModelFits(
    models            = mods,
    dose_levels       = dose_levels,
    posterior         = post,
    simple            = TRUE,
    avg_fit           = TRUE,
    probability_scale = TRUE
  )
  
  # avgFit should be close-ish to p_true at the observed doses when signal is strong
  avg_pred <- fits$avgFit$pred_values
  expect_true(all(avg_pred >= 0 & avg_pred <= 1))
  
  # Probability-scale tolerance: 0.08 is fairly strict with n=250
  expect_equal(as.numeric(avg_pred), as.numeric(p_true), tolerance = 0.08)
  
  # max_effect should be close to true range
  expect_equal(fits$avgFit$max_effect, max(p_true) - min(p_true), tolerance = 0.08)
})


test_that("Bootstrap quantiles: ordered, bounded, and median sits between 20% and 80%", {
  skip_if_not_installed("RBesT")
  skip_if_not_installed("DoseFinding")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  
  set.seed(104)
  
  dose_levels <- c(0, 1, 2, 4)
  p_true      <- c(0.10, 0.20, 0.33, 0.55)
  
  mu_hat <- qlogis(p_true)
  se_hat <- rep(0.35, length(dose_levels))
  
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = mu_hat[i], s = 1.5), sigma = 2)
    }),
    c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
  )
  
  post <- getPosterior(
    prior_list        = prior_list,
    mu_hat            = mu_hat,
    S_hat             = diag(se_hat^2),
    probability_scale = TRUE
  )
  
  mods <- DoseFinding::Mods(
    linear = NULL,
    emax   = c(0.5, 1.2),
    doses  = dose_levels,
    maxEff = 2
  )
  
  fits <- getModelFits(
    models            = mods,
    dose_levels       = dose_levels,
    posterior         = post,
    simple            = TRUE,
    avg_fit           = TRUE,
    probability_scale = TRUE
  )
  
  bq <- getBootstrapQuantiles(
    model_fits        = fits,
    quantiles         = c(0.2, 0.5, 0.8),
    n_samples         = 120,
    doses             = dose_levels,
    probability_scale = TRUE
  )
  
  # Bounds
  expect_true(all(bq$q_val[bq$sample_type == "abs"]  >= 0 & bq$q_val[bq$sample_type == "abs"] <= 1))
  expect_true(all(bq$q_val[bq$sample_type == "diff"] >= -1 & bq$q_val[bq$sample_type == "diff"] <= 1))
  
  # Quantile ordering check (within each model/dose/sample_type)
  # q(0.2) <= q(0.5) <= q(0.8)
  suppressWarnings({
    split_keys <- interaction(bq$model, bq$dose, bq$sample_type, drop = TRUE)
  })
  parts <- split(bq, split_keys)
  for (part in parts) {
    qs <- part$q_prob
    vals <- part$q_val
    o <- order(qs)
    expect_true(all(diff(vals[o]) >= -1e-10))
  }
  
  # Median between 20% and 80% everywhere
  for (part in parts) {
    q20 <- part$q_val[part$q_prob == 0.2]
    q50 <- part$q_val[part$q_prob == 0.5]
    q80 <- part$q_val[part$q_prob == 0.8]
    expect_true(all(q20 <= q50 + 1e-10))
    expect_true(all(q50 <= q80 + 1e-10))
  }
})


test_that("MED: coherent behavior when delta is above/below achievable effect (probability scale)", {
  skip_if_not_installed("RBesT")
  skip_if_not_installed("DoseFinding")
  
  set.seed(105)
  
  dose_levels <- c(0, 1, 2, 4)
  p_true      <- c(0.10, 0.20, 0.35, 0.60)
  
  mu_hat <- qlogis(p_true)
  se_hat <- rep(0.30, length(dose_levels))
  
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = mu_hat[i], s = 1.2), sigma = 2)
    }),
    c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
  )
  
  post <- getPosterior(prior_list = prior_list, mu_hat = mu_hat, S_hat = diag(se_hat^2), probability_scale = TRUE)
  
  mods <- DoseFinding::Mods(
    linear = NULL,
    emax   = c(0.5, 1.2),
    doses  = dose_levels,
    maxEff = 2
  )
  
  fits <- getModelFits(
    models            = mods,
    dose_levels       = dose_levels,
    posterior         = post,
    simple            = TRUE,
    avg_fit           = TRUE,
    probability_scale = TRUE
  )
  
  # Achievable effect from avgFit in probability space
  eff <- fits$avgFit$max_effect
  expect_true(eff >= 0 && eff <= 1)
  
  # delta > effect => should not be reached (for avgFit column at least)
  med_hi <- getMED(delta = min(0.99, eff + 0.15), model_fits = fits, probability_scale = TRUE)
  expect_equal(med_hi["med_reached", "avgFit"], 0)
  
  # delta < effect => should be reached
  med_lo <- getMED(delta = max(0.001, eff - 0.15), model_fits = fits, probability_scale = TRUE)
  expect_equal(med_lo["med_reached", "avgFit"], 1)
})


test_that("Edge case: rare events (0 or all events in a group) yields finite posterior summaries", {
  skip_if_not_installed("RBesT")
  
  set.seed(106)
  
  dose_levels <- c(0, 1, 2)
  n_per_dose  <- 120
  
  # Construct rare-event scenario: control has 0 events, middle has a few, high has many
  y0 <- rep(0, n_per_dose)                              # 0%
  y1 <- rbinom(n_per_dose, 1, 0.03)                     # ~3%
  y2 <- rbinom(n_per_dose, 1, 0.70)                     # high
  
  dat <- data.frame(
    simulation = 1L,
    dose       = rep(dose_levels, each = n_per_dose),
    response   = c(y0, y1, y2)
  )
  attr(dat, "probability_scale") <- TRUE
  
  # Mild prior regularization helps separation-ish scenarios remain finite
  p_center <- c(0.02, 0.05, 0.60)
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = qlogis(p_center[i]), s = 2.0), sigma = 2)
    }),
    c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
  )
  
  post <- getPosterior(prior_list = prior_list, data = dat, probability_scale = TRUE)
  smry <- summary(post)
  
  expect_true(all(is.finite(smry[, "mean"])))
  expect_true(all(is.finite(smry[, "sd"])))
  
  # Posterior mean ordering should reflect data direction (not necessarily strict, but likely)
  # Control logit mean should be smallest; high dose should be largest
  expect_true(smry["Ctr", "mean"] < smry["DG_2", "mean"])
})
