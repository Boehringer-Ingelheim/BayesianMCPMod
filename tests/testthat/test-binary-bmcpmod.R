test_that("performBayesianMCPMod: binary end-to-end returns MCP + fitted models (simple, fast)", {
  skip_if_not_installed("RBesT")
  skip_if_not_installed("DoseFinding")
  
  set.seed(14)
  
  dose_levels <- c(0, 1, 2, 4)
  probs <- c(0.10, 0.20, 0.30, 0.55)
  
  mu_hat <- qlogis(probs)
  se_hat <- rep(0.40, length(dose_levels))
  
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = mu_hat[i], s = 1.4), sigma = 2)
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
  
  n_patients <- c(50, 50, 50, 50)
  
  contr <- getContr(
    mods         = mods,
    dose_levels  = dose_levels,
    dose_weights = n_patients,
    prior_list   = prior_list
  )
  
  # NOTE: your exported entry point is performBayesianMCPMod()
  res <- performBayesianMCPMod(
    posterior_list    = post,
    contr             = contr,
    crit_prob_adj     = 0.5,
    simple            = TRUE,
    avg_fit           = TRUE,
    delta             = 0.10,
    med_selection     = "avgFit",
    n_samples         = 50,
    probability_scale = TRUE
  )
  
  expect_true(is.list(res))
  expect_true(all(c("BayesianMCP", "Mod") %in% names(res)))
  expect_true(is.matrix(res$BayesianMCP) || is.data.frame(res$BayesianMCP))
  expect_true(is.list(res$Mod))
  
  # If delta provided, your function attaches MED as attribute
  expect_true(!is.null(attr(res, "MED")))
})
