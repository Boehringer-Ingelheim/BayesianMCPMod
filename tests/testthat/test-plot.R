test_that("Test plot.modelFits with different model_fits input", {
  library(BayesianMCPMod)
  library(clinDR)
  library(dplyr)
  data("metaData")
  testdata <- as.data.frame(metaData)
  dataset <- filter(testdata, bname == "VIAGRA")
  histcontrol <- filter(dataset, dose == 0, primtime == 12, indication == "ERECTILE DYSFUNCTION")
  
  ## Create MAP Prior
  hist_data <- data.frame(
    trial = c(1, 2, 3, 4),
    est = histcontrol$rslt,
    se = histcontrol$se,
    sd = histcontrol$sd,
    n = histcontrol$sampsize
  )
  
  sd_tot <- with(hist_data, sum(sd * n) / sum(n))
  
  
  gmap <- RBesT::gMAP(
    formula    = cbind(est, se) ~ 1 | trial,
    weights    = hist_data$n,
    data       = hist_data,
    family     = gaussian,
    beta.prior = cbind(0, 100 * sd_tot),
    tau.dist   = "HalfNormal",
    tau.prior  = cbind(0, sd_tot / 4)
  )
  
  prior_ctr <- RBesT::robustify(
    priormix = RBesT::automixfit(gmap),
    weight   = 0.5,
    mean     = 1.4,
    sigma    = sd_tot
  )
  
  # RBesT::ess(prior_ctr)
  
  ## derive prior for treatment
  ## weak informative using same parameters as for robustify component
  prior_trt <- RBesT::mixnorm(
    comp1 = c(w = 1, m = 1.4, n = 1),
    sigma = sd_tot,
    param = "mn"
  )
  dose_levels <- c(0, 50, 100, 200)
  ## combine priors in list
  prior_list <- c(list(prior_ctr), rep(list(prior_trt), times = length(dose_levels[-1])))
  
  # Pre-Specification (B)MCPMod
  
  ## candidate models for MCPMod
  # linear function - no guestimates needed
  exp <- DoseFinding::guesst(
    d = 50,
    p = c(0.2),
    model = "exponential",
    Maxd = max(dose_levels)
  )
  emax <- DoseFinding::guesst(
    d = 100,
    p = c(0.9),
    model = "emax"
  )
  
  
  mods <- DoseFinding::Mods(
    linear      = NULL,
    emax        = emax,
    exponential = exp,
    doses       = dose_levels,
    maxEff      = 10,
    placEff     = 1.4
  )
  
  # Simulation of new trial
  ## Note: This part will be simplified and direct results from one trial will be used
  mods_sim <- DoseFinding::Mods(
    emax        = emax,
    doses       = dose_levels,
    maxEff      = 12,
    placEff     = 1.4
  )
  
  n_patients <- c(50, 50, 50, 50)
  data <- simulateData(
    n_patients  = n_patients,
    dose_levels = dose_levels,
    sd          = sd_tot,
    mods        = mods_sim,
    n_sim       = 1
  )
  
  data_emax <- data[, c("simulation", "dose", "emax")]
  names(data_emax)[3] <- "response"
  
  posterior_emax <- getPosterior(
    data       = data_emax,
    prior_list = prior_list
  )
  
  # Evaluation of Bayesian MCPMod
  
  contr_mat <- DoseFinding::optContr(
    models = mods,
    doses  = dose_levels,
    w      = n_patients
  )
  ## Calculation of critical value can be done with critVal function
  crit_val_equal <- DoseFinding:::critVal(contr_mat$corMat, alpha = 0.05, df = 0, alternative = "one.sided")
  crit_pval <- stats::pnorm(crit_val_equal)
  
  ess_prior <- round(unlist(lapply(prior_list, RBesT::ess)))
  
  ### Evaluation of Bayesian MCPMod
  contr_mat_prior <- DoseFinding::optContr(
    models = mods,
    doses  = dose_levels,
    w      = n_patients + ess_prior
  )
  
  BMCP_result <- performBayesianMCP(
    posterior_list = list(posterior_emax),
    contr          = contr_mat_prior,
    crit_prob      = crit_pval
  )
  
  # Model fit
  # This part is currently not working
  post_observed <- posterior_emax
  model_shapes <- c("linear", "emax", "exponential")
  
  
  # Option a) Simplified approach by using frequentist idea
  fit_simple <- getModelFits(
    models      = model_shapes,
    dose_levels = dose_levels,
    posterior   = post_observed,
    simple      = TRUE
  )
  
  # Option b) Making use of the complete posterior distribution
  fit <- getModelFits(
    models      = model_shapes,
    dose_levels = dose_levels,
    posterior   = post_observed,
    simple      = FALSE
  )
  
  # Test with default parameters and more models
  plot1 <- plot.modelFits(fit)
  expect_s3_class(plot1, "ggplot")
  
  # Test with cr_intv = TRUE and more models
  plot2 <- plot.modelFits(fit, cr_intv = TRUE)
  expect_s3_class(plot2, "ggplot")
  
  # Test with gAIC = FALSE and more models
  plot3 <- plot.modelFits(fit, gAIC = FALSE)
  expect_s3_class(plot3, "ggplot")
  
  # Test with avg_fit = FALSE and more models
  plot4 <- plot.modelFits(fit, avg_fit = FALSE)
  expect_s3_class(plot4, "ggplot")
  
  # Test with all non-default parameters and more models
  plot5 <- plot.modelFits(fit, cr_intv = TRUE, gAIC = FALSE, avg_fit = FALSE)
  expect_s3_class(plot5, "ggplot")
})