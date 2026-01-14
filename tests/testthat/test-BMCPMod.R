# Tests for getCritProb ---------------------------------------------------



# getCritProb relies on DoseFinding, which we assumes works correctly, so the tests here are minimal

test_that("getCritProb returns the right type of value under normal case", {
  crit_pval <- getCritProb(
    mods = mods,
    dose_levels = dose_levels,
    dose_weights = n_patients,
    alpha_crit_val = alpha_crit_val
  )

  expect_type(
    crit_pval, "double"
  )

  expect_true(
    crit_pval >= 0 & crit_pval <= 1
  )
})


# Tests for getContrMat ---------------------------------------------------



# getContrMat relies on DoseFinding, which we assumes works correctly, so the tests here are minimal

test_that("getContrMat returns the right type of object under normal case", {
  contr_mat <- getContr(
    mods = mods,
    dose_levels = dose_levels,
    dose_weights = n_patients,
    prior_list = prior_list
  )

  expect_s3_class(
    contr_mat, "optContr"
  )
})

test_that("getContrMat works as expected", {
  cov_posterior <- diag(sd^2)

  contr_mat_post_sd <- getContr(
    mods          = mods,
    dose_levels   = dose_levels,
    cov_posterior = cov_posterior
  )

  se_new_trial <- c(0.3, 0.7, 0.9, 2.1)
  se_new_trial <- se_new_trial[1:2]

  contr_mat_se_new <- getContr(
    mods = mods,
    dose_levels = dose_levels,
    cov_new_trial = diag(se_new_trial^2)
  )

  # Length mismatch for se_new_trial should error
  expect_error(
    getContr(
      mods = mods,
      dose_levels = dose_levels,
      se_new_trial = se_new_trial[-1]
    )
  )


  expect_s3_class(
    contr_mat_post_sd, "optContr"
  )

  expect_no_error(contr_mat_post_sd)

  expect_s3_class(
    contr_mat_se_new, "optContr"
  )

  expect_no_error(contr_mat_se_new)


  expect_error(
    getContr(
      mods = mods,
      dose_levels = dose_levels
    )
  )
})


# Tests for performBayesianMCP --------------------------------------------



test_that("performBayesianMCP returns the right type of object under normal case", {
  data <- simulateData(
    n_patients  = n_patients,
    dose_levels = dose_levels,
    sd          = sd,
    mods        = mods,
    n_sim       = n_sim
  )

  posterior_list <- getPosterior(
    data = getModelData(data, names(mods)[1]),
    prior_list = prior_list
  )

  contr_mat <- getContr(
    mods = mods,
    dose_levels = dose_levels,
    dose_weights = n_patients,
    prior_list = prior_list
  )

  crit_pval <- getCritProb(
    mods = mods,
    dose_levels = dose_levels,
    dose_weights = n_patients,
    alpha_crit_val = alpha_crit_val
  )

  b_mcp <- performBayesianMCP(
    posterior_list = posterior_list,
    contr = contr_mat,
    crit_prob_adj = crit_pval
  )

  expect_s3_class(
    b_mcp,
    "BayesianMCP"
  )

  expect_true(
    attr(b_mcp, "critProbAdj") == crit_pval
  )

  expect_type(
    attr(b_mcp, "essAvg"), "logical"
  )

  expect_type(
    attr(b_mcp, "successRate"), "double"
  )


  expect_type(b_mcp, "double")
})


# Tests for performBayesianMCPMod -----------------------------------------



test_that("performBayesianMCPMod returns the right type of object under normal case", {
  b_mcp_mod <- performBayesianMCPMod(
    posterior_list = posterior_list,
    contr = contr_mat,
    crit_prob_adj = crit_pval
  )

  expect_s3_class(
    b_mcp_mod,
    "BayesianMCPMod"
  )

  expect_true(
    all(names(b_mcp_mod) == c("BayesianMCP", "Mod"))
  )
})


# Tests for addSignificance -----------------------------------------------


test_that("addSignificance attaches flags per model and validates input length", {
  addSignificance_fn <- tryCatch(
    getFromNamespace("addSignificance", "BayesianMCPMod"),
    error = function(e) NULL
  )
  skip_if(is.null(addSignificance_fn), "addSignificance not exported/available")

  models <- c("emax", "linear")
  dose_levels <- c(0, 1, 2, 4, 8)

  # Strongly convex pattern: tiny effects at low/mid doses, big jump at the top dose.
  posterior_list <- list(
    Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0.0, s = 0.7),  sigma = 1.2),
    DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 0.6, s = 0.7),  sigma = 1.2),
    DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.4, s = 0.7),  sigma = 1.2),
    DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 3.4, s = 0.7),  sigma = 1.2),
    DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 7.8, s = 0.7),  sigma = 1.2)
  )


  fit <- getModelFits(
    models     = models,
    posterior  = posterior_list,
    dose_levels = dose_levels
  )


  # Flags length matches -> flags should be attached per entry
  out <- addSignificance_fn(fit, c(TRUE, FALSE))
  expect_true(is.list(out) && all(names(out) == names(fit)))
  expect_false(out$linear$significant.linear)
  expect_false(out$emax$significant.emax)

  # Mismatched length should raise an error
  expect_error(addSignificance_fn(fit, list(TRUE)))
})




