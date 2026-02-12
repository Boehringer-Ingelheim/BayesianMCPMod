test_that("simulateData: binary mode produces 0/1 outcomes and stores attribute", {
  skip_if_not_installed("DoseFinding")
  skip_if_not_installed("RBesT")
  
  set.seed(1)
  
  dose_levels <- c(0, 1, 2, 4)
  n_patients  <- c(30, 30, 30, 30)
  
  mods <- DoseFinding::Mods(
    linear = NULL,
    emax   = c(0.5, 1.2),
    doses  = dose_levels,
    maxEff = 2
  )
  
  sim <- simulateData(
    n_patients        = n_patients,
    dose_levels       = dose_levels,
    sd                = NULL,
    mods              = mods,
    n_sim             = 3,
    probability_scale = TRUE
  )
  
  expect_true(is.data.frame(sim))
  expect_true(all(c("simulation", "ptno", "dose") %in% names(sim)))
  expect_true(attr(sim, "probability_scale"))
  
  # All response columns should be 0/1
  resp_cols <- setdiff(names(sim), c("simulation", "ptno", "dose"))
  expect_true(length(resp_cols) >= 1)
  
  for (cn in resp_cols) {
    expect_true(all(sim[[cn]] %in% c(0, 1)))
  }
  
  # Must error if continuous mode but sd is missing
  expect_error(
    simulateData(
      n_patients        = n_patients,
      dose_levels       = dose_levels,
      sd                = NULL,
      mods              = mods,
      n_sim             = 1,
      probability_scale = FALSE
    ),
    regexp = "Must provide 'sd'|sd"
  )
})

test_that("simulateData: binary mode works with custom dr_means and sd=NULL", {
  skip_if_not_installed("RBesT")
  
  set.seed(2)
  
  dose_levels <- c(0, 1, 2, 4)
  n_patients  <- c(25, 25, 25, 25)
  
  # dr_means are on the *logit* scale (inv_logit used internally)
  dr_means_logit <- c(-2, -1, 0, 1)
  
  sim <- simulateData(
    n_patients        = n_patients,
    dose_levels       = dose_levels,
    sd                = NULL,
    mods              = NULL,
    n_sim             = 2,
    dr_means          = dr_means_logit,
    probability_scale = TRUE
  )
  
  resp_cols <- setdiff(names(sim), c("simulation", "ptno", "dose"))
  expect_equal(length(resp_cols), 1L)
  expect_true(all(sim[[resp_cols]] %in% c(0, 1)))
})
