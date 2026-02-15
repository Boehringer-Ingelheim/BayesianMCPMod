test_that("getModelFits: probability_scale=TRUE gives pred_values in [0,1] and predict() is consistent", {
  skip_if_not_installed("RBesT")
  skip_if_not_installed("DoseFinding")
  
  set.seed(11)
  
  dose_levels <- c(0, 1, 2, 4)
  n_per_dose  <- 60
  probs       <- c(0.10, 0.20, 0.30, 0.55)
  
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
  
  expect_s3_class(fits, "modelFits")
  expect_true(isTRUE(attr(fits, "probability_scale")))
  
  for (nm in names(fits)) {
    pv <- fits[[nm]]$pred_values
    expect_true(is.numeric(pv))
    expect_true(all(pv >= 0 & pv <= 1))
    expect_equal(fits[[nm]]$max_effect, max(pv) - min(pv), tolerance = 1e-12)
  }
  
  # predict() should match inv_logit of logit-scale predictions
  doses_new <- sort(unique(c(0, 0.5, 1, 3, 4)))
  pred_prob <- predict(fits, doses = doses_new, probability_scale = TRUE)
  pred_log  <- predict(fits, doses = doses_new, probability_scale = FALSE)
  
  for (nm in names(pred_prob)) {
    expect_equal(
      as.numeric(pred_prob[[nm]]),
      as.numeric(RBesT::inv_logit(pred_log[[nm]])),
      tolerance = 1e-10
    )
  }
})
